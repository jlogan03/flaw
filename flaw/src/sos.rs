use num_traits::{MulAdd, Num};

use crate::{AlignedArray, Ring};

/// Single-Input-Single-Output, cascaded Second Order Sections filter.
#[derive(Clone, Copy)]
pub struct SisoSosFilter<const ORDER: usize, T: Num + Copy + MulAdd<Output = T>> {
    /// Latest output
    y: T,
    /// Internal state: two filter delays for each section
    z: [[T; 2]; ORDER],
    /// SOS coefficients ordered as [b0, b1, b2, a1, a2] per section.
    /// These correspond to a filter transfer function of:
    ///    H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
    /// Note that some SOS implementations (e.g. SciPy) use six coefficients,
    /// including a0, but here a0 is assumed to be unity, and is omitted to reduce memory usage.
    sos: [[T; 5]; ORDER],
}

impl<const ORDER: usize, T: Num + Copy + MulAdd<Output = T>> SisoSosFilter<ORDER, T> {

    #[inline]
    pub fn update(&mut self, u: T) -> T {
        let mut input = u; // Input to each section
        let mut output = T::zero(); // Output of each section
        for s in 0..ORDER {
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
        let mut sos_ = [[T::zero(); 5]; ORDER];
        sos_.copy_from_slice(sos);

        Self {
            y: T::zero(),
            z: [[T::zero(); 2]; ORDER],
            sos: sos_,
        }
    }
}
