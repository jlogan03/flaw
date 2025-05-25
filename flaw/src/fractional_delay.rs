//! Fractional delay all-pass filter, typically used for aligning scanned samples to a
//! single effective sample time.
//!
//! ## References
//!  
//! \[1\] “Lagrange Interpolation | Physical Audio Signal Processing.” Accessed: May 14, 2025. [Online]. Available: https://www.dsprelated.com/freebooks/pasp/Lagrange_Interpolation.html
//!
//! \[2\] “Lagrange polynomial,” Wikipedia. Apr. 16, 2025. Accessed: May 14, 2025. [Online]. Available: https://en.wikipedia.org/wiki/Lagrange_polynomial

use crate::SisoFirFilter;
use num_traits::Num;


/// Calculate taps for a Lagrange polynomial fractional delay filter
/// which creates a linear-phase lag of `delay` as a fraction of one sample period.
/// For best results, `delay` should be between 0 and ORDER - 1; typical applications
/// use `delay` between 0 and 1.
/// 
/// Taps are ordered from most recent sample to least recent sample.
/// 
/// 2 <= ORDER <= 255 is required. The actual maximum usable order varies by
/// numerical field type, and is typically around 10.
pub fn lagrange_fractional_delay_taps<const ORDER: usize, T: Num + From<u8> + Copy>(
    delay: T,
) -> [T; ORDER] {
    let mut taps = [T::zero(); ORDER];

    const {
        assert!(
            ORDER <= u8::MAX as usize,
            "Filter order must be less than u8::MAX"
        );
        assert!(ORDER >= 2, "Filter order less than 2 is meaningless");
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

    // Enforce that taps sum to exactly 1.0
    let tapsum = taps.iter().fold(T::zero(), |acc, &x| acc + x);
    taps.iter_mut().for_each(|x| *x = *x / tapsum);

    taps
}


/// Generate a fractional-delay filter using Lagrange polynomials, which provides a controlled
/// time delay as a fraction of one sample period, with a _filter_ order of `ORDER`
/// and _polynomial_ order of `ORDER-1`. Maximum order for this function is 255, which is
/// already much larger than practical. Typical order is <10.
///
/// This is an FIR filter whose taps are a baked result of evaluating a minimum-order
/// polynomial at a specific location - just a way of reducing a polynomial interpolation
/// to a dot product, similar to how finite difference stencil coefficients are generated,
/// but in this special case with equal sample spacing, we do not need to invert a matrix.
///
/// Building the taps requires 4 * N^2 - N operations for N taps.
pub fn lagrange_fractional_delay_filter<const ORDER: usize, T: Num + From<u8> + Copy>(
    delay: T,
) -> SisoFirFilter<ORDER, T> {
    let taps: [T; ORDER] = lagrange_fractional_delay_taps(delay);
    SisoFirFilter::new(&taps)
}


#[cfg(test)]
mod test {
    use super::lagrange_fractional_delay_filter;

    /// Check that the fractional delay produces the correct value for special cases
    /// where the signal is a polynomial that can be reproduced exactly by the filter.
    #[test]
    fn test_lagrange_fractional_delay_filter() {

        // Linear signal
        let vals = [4.0, 3.0, 2.0, 1.0];
        //   Linear filter
        let mut filter: crate::SisoFirFilter<2, f64> = lagrange_fractional_delay_filter(0.5);
        vals.iter().for_each(|v| {filter.update(*v);});
        assert_eq!(filter.y(), 1.5);
        //   Quadratic filter
        let mut filter: crate::SisoFirFilter<3, f64> = lagrange_fractional_delay_filter(0.5);
        vals.iter().for_each(|v| {filter.update(*v);});
        assert_eq!(filter.y(), 1.5);
        //   Cubic filter
        let mut filter: crate::SisoFirFilter<4, f64> = lagrange_fractional_delay_filter(0.5);
        vals.iter().for_each(|v| {filter.update(*v);});
        assert_eq!(filter.y(), 1.5);

        // Cubic signal
        let func = |x: f64| {0.3 + 0.18 * x - 0.5 * x.powf(2.0) + 0.7 * x.powf(3.0)};
        let vals = [func(0.0), func(1.0), func(2.0), func(3.0)];
        let mut filter: crate::SisoFirFilter<4, f64> = lagrange_fractional_delay_filter(0.2);
        vals.iter().for_each(|v| {filter.update(*v);});
        let err = (filter.y() - func(2.8)).abs();
        assert!(err < 1e-14);
    }
}