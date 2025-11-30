//! Second Order Sections (SOS) filters.

// Required for float conversions from 64 to 32 bit and log10 on no-std targets
#[cfg(not(feature = "std"))]
use num_traits::Float;

use core::ops::Neg;
use crunchy::unroll;
use num_traits::{FromPrimitive, MulAdd, Num, ToPrimitive};

pub mod tables;
pub use tables::butter2::butter2;
pub use tables::butter4::butter4;
pub use tables::butter6::butter6;

use crate::AlignedArray;

#[cfg(test)]
mod test_helpers;

/// Single-Input-Single-Output, cascaded Second Order Sections filter.
#[derive(Clone, Copy)]
pub struct SisoSosFilter<
    const SECTIONS: usize,
    T: Num + Copy + MulAdd<Output = T> + Neg<Output = T>,
> {
    /// Latest output
    y: T,
    /// Internal state: two filter delays for each section
    z: AlignedArray<[T; 2], SECTIONS>,
    /// SOS coefficients ordered as [b0, b1, b2, a1, a2] per section.
    /// These correspond to a filter transfer function of:
    ///
    ///    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
    ///
    /// Note that some SOS implementations (e.g. SciPy) use six coefficients,
    /// including a0, but here a0 is assumed to be unity, and is omitted to reduce memory usage.
    sos: AlignedArray<[T; 5], SECTIONS>, // Using AlignedArray did not measurably change performance on an i7-8550U CPU, but might help on other platforms
}

impl<const SECTIONS: usize, T> SisoSosFilter<SECTIONS, T>
where
    T: Num + Copy + MulAdd<Output = T> + Neg<Output = T> + FromPrimitive + ToPrimitive, // FromPrimitive needed for conversion from f64 in SOS coef lookup tables.
{
    /// Evaluate the next estimated value based on the latest measurement
    /// in 9N floating-point ops for a filter with N sections (order up to 2*N).
    #[inline]
    pub fn update(&mut self, u: T) -> T {
        let mut input = u; // Input to each section
        let mut output = T::zero(); // Output of each section

        // `crunchy::unroll` requires a literal upper bound; use a conservative
        // limit and let the assert catch unsupported configurations.
        const {
            assert!(
                SECTIONS <= 8,
                "SOS filters with SECTIONS > 8 are not supported"
            );
        }

        // Unrolling this loop gives a ~5% speedup of the sos butter4 f64 benchmark
        // on a i7-8550U CPU.
        unroll! {
            for s < 8 in 0..SECTIONS {
                let b0 = self.sos[s][0];
                let b1 = self.sos[s][1];
                let b2 = self.sos[s][2];
                let a1 = self.sos[s][3];
                let a2 = self.sos[s][4];

                #[cfg(not(feature = "fma"))]
                {
                    // Direct Form II Transposed implementation
                    output = b0 * input + self.z[s][0];
                    // Update the filter delays, they will be used the next time this function is called
                    self.z[s][0] = b1 * input - a1 * output + self.z[s][1];
                    self.z[s][1] = b2 * input - a2 * output;
                }

                // The FMA implementation is ~40% faster, with target-cpu=x86-64-v3 on a i7-8550U CPU.
                #[cfg(feature = "fma")]
                {
                    // Direct Form II Transposed implementation
                    output = b0.mul_add(input, self.z[s][0]);  // b0 * input + self.z[s][0]
                    // Update the filter delays, they will be used the next time this function is called
                    self.z[s][0] = b1.mul_add(input, a1.mul_add(-output, self.z[s][1])); // b1 * input - a1 * output + self.z[s][1]
                    self.z[s][1] = b2.mul_add(input, -a2 * output); // b2 * input - a2 * output
                }

                // Cascaded sections: output of this section is input to next
                #[allow(unused_assignments)] // rustc warns because the assignment to input is unused on the last section
                {
                    input = output;
                }
            }
        }
        // Overall output of the filter is the output of the last sections
        self.y = output;
        self.y
    }

    /// Reset internal state to zero.
    pub fn reset(&mut self) {
        self.y = T::zero();
        self.z = AlignedArray([[T::zero(); 2]; SECTIONS]);
    }

    /// Set filter internal state to the steady value
    /// achieved for input `u`. For filters with unity steady-state gain,
    /// this will also produce an output reading of `u`.
    pub fn set_steady_state(&mut self, u: T) -> Result<(), &'static str> {
        let mut input = u;
        let mut overall_ss_gain = 1.0;
        for s in 0..SECTIONS {
            let b0 = self.sos[s][0];
            let b1 = self.sos[s][1];
            let b2 = self.sos[s][2];
            let a1 = self.sos[s][3];
            let a2 = self.sos[s][4];

            // Calculate the steady-state output of this section
            let input_f64 = input.to_f64().ok_or("Conversion to f64 failed")?;
            let section_f64: [f64; 5] = [
                b0.to_f64().ok_or("Conversion to f64 failed")?,
                b1.to_f64().ok_or("Conversion to f64 failed")?,
                b2.to_f64().ok_or("Conversion to f64 failed")?,
                a1.to_f64().ok_or("Conversion to f64 failed")?,
                a2.to_f64().ok_or("Conversion to f64 failed")?,
            ];
            let ss_gain = steady_state_gain_sos(&section_f64);
            overall_ss_gain *= ss_gain; // accumulate the gain of each section
            let output: T = T::from_f64(input_f64 * ss_gain).ok_or("Conversion from f64 failed")?;

            // Set the internal states based on the steady state input and output of this section
            self.z[s][1] = b2 * input - a2 * output;
            self.z[s][0] = b1 * input - a1 * output + self.z[s][1];

            // Cascaded sections: output of this section is input to next
            input = output;
        }

        if (overall_ss_gain - 1.0).abs() < 1e-6 {
            // Try updating the filter and verify that the output matches the input
            debug_assert!(
                // use debug_assert so that release builds can be panic-free
                (self.update(u) - u)
                    .to_f64()
                    .ok_or("Conversion to f64 failed")?
                    .abs()
                    < 1e-6,
            );
        }

        Ok(())
    }

    pub fn new(sos: &[[T; 5]]) -> Self {
        let mut sos_ = [[T::zero(); 5]; SECTIONS];
        sos_.copy_from_slice(sos);

        Self {
            y: T::zero(),
            z: AlignedArray([[T::zero(); 2]; SECTIONS]),
            sos: AlignedArray(sos_),
        }
    }

    /// Build a new low-pass with coefficients interpolated on baked tables.
    pub fn new_interpolated(
        cutoff_ratio: f64,
        log10_cutoff_ratio_grid: &[f64],
        sos_tables: [[&[f64]; 5]; SECTIONS],
    ) -> Result<Self, &'static str> {
        let log10_cutoff_ratio = cutoff_ratio.log10();

        // Check table bounds
        let mut extrapolated = [false; 1];
        interpn::multicubic::rectilinear::check_bounds(
            &[log10_cutoff_ratio_grid],
            &[&[log10_cutoff_ratio]],
            1e-6,
            &mut extrapolated,
        )?;

        if extrapolated[0] {
            return Err("Selected cutoff ratio is outside the grid");
        }

        let mut sos = [[T::zero(); 5]; SECTIONS];
        for sec in 0..SECTIONS {
            for coeff in 0..5 {
                let val_f64: f64 = interpn::MulticubicRectilinear::<'_, _, 1>::new(
                    &[log10_cutoff_ratio_grid],
                    sos_tables[sec][coeff],
                    true,
                )?
                .interp_one([log10_cutoff_ratio])?;
                sos[sec][coeff] = T::from_f64(val_f64).ok_or("Conversion from f64 failed")?;
            }
        }

        // Correct the DC gain of the filter to 1.
        // First, calculate the DC gain we'd get with the coefficients as is after interpolation.
        let mut dc_gain: f64 = 1.0;
        for section in sos.iter() {
            let b0 = section[0].to_f64().ok_or("Conversion to f64 failed")?;
            let b1 = section[1].to_f64().ok_or("Conversion to f64 failed")?;
            let b2 = section[2].to_f64().ok_or("Conversion to f64 failed")?;
            let a1 = section[3].to_f64().ok_or("Conversion to f64 failed")?;
            let a2 = section[4].to_f64().ok_or("Conversion to f64 failed")?;
            let sec_dc_gain = (b0 + b1 + b2) / (1.0 + a1 + a2);
            dc_gain *= sec_dc_gain;
        }
        // Now scale the numerator coefficients to get unity gain at DC.
        // Apply the required scaling in even parts to each section.
        let correction = T::from_f64(dc_gain.powf(-1.0 / SECTIONS as f64))
            .ok_or("Conversion from f64 failed")?;
        for section in sos.iter_mut() {
            section[0] = section[0] * correction;
            section[1] = section[1] * correction;
            section[2] = section[2] * correction;
        }
        // Verify the DC gain after scaling.
        dc_gain = 1.0;
        for section in sos.iter() {
            let b0 = section[0].to_f64().ok_or("Conversion to f64 failed")?;
            let b1 = section[1].to_f64().ok_or("Conversion to f64 failed")?;
            let b2 = section[2].to_f64().ok_or("Conversion to f64 failed")?;
            let a1 = section[3].to_f64().ok_or("Conversion to f64 failed")?;
            let a2 = section[4].to_f64().ok_or("Conversion to f64 failed")?;
            let sec_dc_gain = (b0 + b1 + b2) / (1.0 + a1 + a2);
            dc_gain *= sec_dc_gain;
        }
        if (dc_gain - 1.0).abs() > 1e-6 {
            return Err("DC gain correction failed");
        }

        Ok(Self::new(&sos))
    }
}

/// Calculate the steady-state gain of a single Second Order Section.
pub fn steady_state_gain_sos(section: &[f64; 5]) -> f64 {
    let b0 = section[0];
    let b1 = section[1];
    let b2 = section[2];
    let a1 = section[3];
    let a2 = section[4];

    // The formula for steady state gain is the z-transform evaluated at z = 1.
    // Derivation: the steady state gain is the final value of the step response.
    // Let H(z) be the z-transform of the second order section.
    // The z-transform of the step response is:
    //   Y(z) = (1 / (1 - z^-1)) * H(z)
    // The final value theorem for z-transforms states:
    //   lim_{n->inf} y[n] = lim_{z->1} (z - 1) * Y(z)
    // Applying this to the step response gives:
    //   lim_{n->inf} y[n] = lim_{z->1} (z - 1) / (1 - z^-1) * H(z)
    // The limit of the first term is 1, so we have:
    //   lim_{n->inf} y[n] = H(1)
    // assuming that H(z) is stable so the limit exists.
    // See https://en.wikipedia.org/wiki/Z-transform#Properties
    (b0 + b1 + b2) / (1.0 + a1 + a2)
}

#[cfg(feature = "std")]
#[cfg(test)]
mod test {
    use super::SisoSosFilter;
    use super::test_helpers::simulate_gain_sinewave;

    #[test]
    fn test_sos_butter_gain() {
        // Coefficients for a 4th order lowpass Butterworth filter with fc/fs = 0.05.
        // Coefficients computed with scipy.signal.butter.
        let mut filter = SisoSosFilter::<2, f64>::new(&[
            [
                4.16599204e-04,
                8.33198409e-04,
                4.16599204e-04,
                -1.47967422e+00,
                5.55821543e-01,
            ],
            [
                1.00000000e+00,
                2.00000000e+00,
                1.00000000e+00,
                -1.70096433e+00,
                7.88499740e-01,
            ],
        ]);
        // Measure the gain of the filter on sine waves of different frequencies.
        // Compare the measured gains to the expected gains. The expected gains were
        // calculated with scipy.signal.freqz_sos in scripts/sos_test.py.
        // f is frequency normalized to the sample frequency.
        for (f, gain_expected) in [
            (0.010, 1.000),                // gain should be 1 well below the cutoff frequency
            (0.050, 1.0 / f64::sqrt(2.0)), // gain should be ~1/sqrt(2) at the cutoff frequency
            (0.100, 5.643e-2),             // gain should be << 1 well above the cutoff frequency
        ] {
            let gain = simulate_gain_sinewave(&mut filter, f, 1024);
            let err = (gain - gain_expected).abs();
            assert!(
                err < 0.01,
                "f = {f}: gain = {gain}, expected {gain_expected}, err = {err}"
            );
        }
    }

    #[test]
    fn test_set_steady_state() {
        // Coefficients for a 4th order lowpass Butterworth filter with fc/fs = 0.05.
        // Coefficients computed with scipy.signal.butter.
        let mut filter = SisoSosFilter::<2, f64>::new(&[
            [
                4.16599204e-04,
                8.33198409e-04,
                4.16599204e-04,
                -1.47967422e+00,
                5.55821543e-01,
            ],
            [
                1.00000000e+00,
                2.00000000e+00,
                1.00000000e+00,
                -1.70096433e+00,
                7.88499740e-01,
            ],
        ]);

        for u in [-1.1, 0.0, 0.5, 1.0, 5.59823] {
            filter.set_steady_state(u).unwrap();
            let y = filter.update(u);
            let err = (y - u).abs();
            assert!(err < 1e-6, "u = {u}: y = {y}, expected {u}, err = {err}");
        }
    }
}
