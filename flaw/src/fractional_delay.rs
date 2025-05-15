//! Fractional delay filter using Lagrange polynomials
//! 
//! ## References
//!  
//! \[1\] “Lagrange Interpolation | Physical Audio Signal Processing.” Accessed: May 14, 2025. [Online]. Available: https://www.dsprelated.com/freebooks/pasp/Lagrange_Interpolation.html
//! 
//! \[2\] “Lagrange polynomial,” Wikipedia. Apr. 16, 2025. Accessed: May 14, 2025. [Online]. Available: https://en.wikipedia.org/wiki/Lagrange_polynomial

use crate::SisoFirFilter;
use num_traits::Num;

/// Generate a fractional-delay filter using Lagrange polynomials, which provides a controlled
/// time delay as a fraction of one sample period, with a _filter_ order of `ORDER`
/// and _polynomial_ order of `ORDER-1`. Maximum order for this function is 255, which is
/// already larger than practical. Typical order is <10.
///
/// This is an FIR filter whose taps are a baked result of evaluating a minimum-order
/// polynomial at a specific location - just a way of reducing a polynomial interpolation
/// to a dot product, similar to how finite difference stencil coefficients are generated,
/// but in this special case with equal sample spacing, we do not need to invert a matrix.
/// 
/// Building the taps requires 4 * N^2 - N operations for N taps.
pub fn fractional_delay_taps<const ORDER: usize, T: Num + From<u8> + Copy>(delay: T) -> [T; ORDER] {
    let mut taps = [T::zero(); ORDER];

    const {
        assert!(ORDER <= u8::MAX as usize, "Filter order must be less than u8::MAX")
    }

    for k in 0..ORDER as u8 {
        let mut coeff = T::one();
        let kv = T::from(k);
        for m in 0..ORDER as u8 {
            let mv = T::from(m);
            if m != k {
                coeff = coeff * (delay - mv) / (kv - mv);
            }
        }
        taps[k as usize] = coeff;
    }

    // TODO: enforce that taps sum to exactly 1.0

    taps
}
