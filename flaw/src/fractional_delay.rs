use crate::SisoFirFilter;
use num_traits::Float;

/// Generate a fractional-delay aka Farrow filter, which provides a controlled
/// time delay as a fraction of one sample period, with a _filter_ order of ORDER
/// and _polynomial_ order of ORDER-1.
///
/// This is an FIR filter whose taps are a baked result of evaluating a minimum-order
/// polynomial at a specific location - just a way of reducing a polynomial interpolation
/// to a dot product, similar to how finite difference stencil coefficients are generated,
/// but in this special case with equal sample spacing, we do not need to invert a matrix.
pub fn fractional_delay_taps<const ORDER: usize, T: Float + Copy>(delay: T) -> [T; ORDER] {
    let mut taps = [T::zero(); ORDER];

    for k in 0..ORDER {
        let mut coeff = T::one();
        let kv = T::from(k).unwrap();
        for m in 0..ORDER {
            let mv = T::from(m).unwrap();
            if m != k {
                coeff = coeff * (delay - mv) / (kv - mv);
            }
        }
        taps[k] = coeff;
    }

    taps
}
