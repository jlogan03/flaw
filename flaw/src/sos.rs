use num_traits::{FromPrimitive, MulAdd, Num};

/// Single-Input-Single-Output, cascaded Second Order Sections filter.
#[derive(Clone, Copy)]
pub struct SisoSosFilter<const SECTIONS: usize, T: Num + Copy + MulAdd<Output = T>> {
    /// Latest output
    y: T,
    /// Internal state: two filter delays for each section
    z: [[T; 2]; SECTIONS],
    /// SOS coefficients ordered as [b0, b1, b2, a1, a2] per section.
    /// These correspond to a filter transfer function of:
    ///    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
    /// Note that some SOS implementations (e.g. SciPy) use six coefficients,
    /// including a0, but here a0 is assumed to be unity, and is omitted to reduce memory usage.
    sos: [[T; 5]; SECTIONS],
}

impl<const SECTIONS: usize, T: Num + Copy + MulAdd<Output = T>> SisoSosFilter<SECTIONS, T> {
    /// Evaluate the next estimated value based on the latest measurement
    /// in 9N floating-point ops for a filter with N sections (order up to 2*N).
    #[inline]
    pub fn update(&mut self, u: T) -> T {
        let mut input = u; // Input to each section
        let mut output = T::zero(); // Output of each section
        for s in 0..SECTIONS {
            let b0 = self.sos[s][0];
            let b1 = self.sos[s][1];
            let b2 = self.sos[s][2];
            let a1 = self.sos[s][3];
            let a2 = self.sos[s][4];

            // Direct Form II Transposed implementation
            output = b0 * input + self.z[s][0];
            // Update the filter delays, they will be used the next time this function is called
            self.z[s][0] = b1 * input - a1 * output + self.z[s][1];
            self.z[s][1] = b2 * input - a2 * output;
            // Cascaded sections: output of this section is input to next
            input = output;
        }
        // Overall output of the filter is the output of the last sections
        self.y = output;
        self.y
    }

    /// Reset internal state to zero.
    pub fn reset(&mut self) {
        self.y = T::zero();
        self.z = [[T::zero(); 2]; SECTIONS];
    }

    pub fn new(sos: &[[T; 5]]) -> Self {
        let mut sos_ = [[T::zero(); 5]; SECTIONS];
        sos_.copy_from_slice(sos);

        Self {
            y: T::zero(),
            z: [[T::zero(); 2]; SECTIONS],
            sos: sos_,
        }
    }
}

impl<const SECTIONS: usize, T> SisoSosFilter<SECTIONS, T>
where
    T: Num + Copy + MulAdd<Output = T> + FromPrimitive, // FromPrimitive needed for conversion from f64 in SOS coef lookup tables.
{
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
        // TODO is scaling needed here for unity gain at DC?
        Ok(Self::new(&sos))

    }
}
 

#[cfg(feature = "std")]
#[cfg(test)]
 mod test {
    use super::SisoSosFilter;

    #[test]
    fn test_sos_butter_gain() {
        // Coefficients for a 4th order lowpass Butterworth filter with fc/fs = 0.05.
        // Coefficients computed with scipy.signal.butter.
        let mut filter = SisoSosFilter::<2, f64>::new(&[
            [4.16599204e-04, 8.33198409e-04, 4.16599204e-04, -1.47967422e+00, 5.55821543e-01],
            [1.00000000e+00, 2.00000000e+00, 1.00000000e+00, -1.70096433e+00, 7.88499740e-01],
        ]);
        // Measure the gain of the filter on sine waves of different frequencies.
        // Compare the measured gains to the expected gains. The expected gains were
        // calculated with scipy.signal.freqz_sos in scripts/sos_test.py.
        // f is frequency normalized to the sample frequency.
        for (f, gain_expected) in [
            (0.010, 1.000), // gain should be 1 well below the cutoff frequency
            (0.050, 7.057e-1), // gain should be ~1/sqrt(2) at the cutoff frequency
            (0.100, 5.643e-2), // gain should be << 1 well above the cutoff frequency
        ] {
            // Setup
            let n = 1024;
            let input: Vec<f64> = (0..n).map(|i| (2.0 * std::f64::consts::PI * f * i as f64).sin()).collect();
            let original_rms = (input.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
            let mut output = Vec::with_capacity(n);

            // Action
            filter.reset();
            for &u in &input {
                output.push(filter.update(u));
            }

            // Verification
            let output_rms = (output.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
            let gain = output_rms / original_rms;
            let err = (gain - gain_expected).abs();
            assert!(err < 0.01, "f = {f}: gain = {gain}, expected {gain_expected}, err = {err}");
        }
    }
 }