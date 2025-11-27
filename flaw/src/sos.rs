use num_traits::{MulAdd, Num};

/// Single-Input-Single-Output, cascaded Second Order Sections filter.
#[derive(Clone, Copy)]
pub struct SisoSosFilter<const HALF_ORDER: usize, T: Num + Copy + MulAdd<Output = T>> {
    /// Latest output
    y: T,
    /// Internal state: two filter delays for each section
    z: [[T; 2]; HALF_ORDER],
    /// SOS coefficients ordered as [b0, b1, b2, a1, a2] per section.
    /// These correspond to a filter transfer function of:
    ///    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
    /// Note that some SOS implementations (e.g. SciPy) use six coefficients,
    /// including a0, but here a0 is assumed to be unity, and is omitted to reduce memory usage.
    sos: [[T; 5]; HALF_ORDER],
}

impl<const HALF_ORDER: usize, T: Num + Copy + MulAdd<Output = T>> SisoSosFilter<HALF_ORDER, T> {
    /// Evaluate the next estimated value based on the latest measurement
    /// in 9N floating-point ops for a filter with N sections (order 2*N).
    #[inline]
    pub fn update(&mut self, u: T) -> T {
        let mut input = u; // Input to each section
        let mut output = T::zero(); // Output of each section
        for s in 0..HALF_ORDER {
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

    pub fn new(sos: &[[T; 5]]) -> Self {
        let mut sos_ = [[T::zero(); 5]; HALF_ORDER];
        sos_.copy_from_slice(sos);

        Self {
            y: T::zero(),
            z: [[T::zero(); 2]; HALF_ORDER],
            sos: sos_,
        }
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
            let signal: Vec<f64> = (0..n).map(|i| (2.0 * std::f64::consts::PI * f * i as f64).sin()).collect();
            let original_rms = (signal.iter().map(|v| v * v).sum::<f64>() / n as f64).sqrt();
            let mut output = Vec::with_capacity(n);

            // Action
            filter.z = [[0.0; 2]; 2]; // Reset filter state
            for &u in &signal {
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