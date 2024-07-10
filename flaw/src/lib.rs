//! # flaw
//! Control-law-inspired embedded signal filtering, no-std and no-alloc compatible.
//!
//! This library provides a simple method for initializing and updating single-input,
//! single-output infinite-impulse-response filters using 32-bit floats, as well as
//! tabulated filter coefficients for some common filters. Filters evaluate in
//! 4N+1 floating-point operations for a filter of order N.
//!
//! The name `flaw` is short for filter-law, but also refers to the fact that
//! digital IIR filtering with small floating-point types is an inherently flawed
//! approach, in that higher-order and lower-cutoff filters produce very small
//! coefficients that result in floating-point roundoff error. This library makes
//! an attempt to mitigate this problem by providing filter coefficients for a tested
//! domain of validity. The result is a limited, but useful, range of operation
//! where these filters can achieve both accuracy and performance as well
//! as be formulated and initialized in an embedded environment.
//!
//! ## Example: Second-Order Butterworth Filter
//!
//! ```rust
//! // First, choose a cutoff frequency as a fraction of sampling frequency
//! let cutoff_ratio = 1e-3;
//!
//! // Initialize a filter, interpolating coefficients to that cutoff ratio.
//! let mut filter = flaw::butter2(cutoff_ratio).unwrap();  // Errors if extrapolating
//!
//! // Update the filter with a new raw measurement
//! let measurement = 0.3145; // Some number
//! let estimate = filter.update(measurement);  // Latest state estimate
//! ```
//!
//! ## Development Status: Early Days
//!
//! This is in an experimental stage - it appears to work well, but is not fully-validated
//! or fully-featured.
//!
//! * More software testing is needed to guarantee filter performance at interpolated cutoff ratios
//! * More hardware/firmware testing is needed to examine performance on actual microcontrollers
//! * More filter types can be added
//!
//! ## Coefficient Tables
//!
//! Tabulated filters are tested to enforce
//!
//! * <0.1% error in converged step response at the minimum cutoff frequency
//! * <1ppm error in converged step response at the maximum cutoff frequency
//! * <5% error to -3dB attenuation of a sine input at the cutoff frequency at the maximum cutoff ratio
//!   * This error appears to be mainly an issue of discretization in test cases, and could be reduced
//!     by using a better method for testing (fit a sine curve to the result or do gradient-descent
//!     on a cubic interpolator)
//!
//! Each filter with tabulated coefficients has a minimum and maximum cutoff ratio.
//! The minimum value is determined by floating-point error in convergence of a
//! step response, while the maximum value is determined by the accuracy of attenuation
//! at the cutoff frequency as the cutoff ratio approaches the Nyquist frequency.
//!
//! Coefficients for a given filter are interpolated on these tables using a
//! cubic Hermite method with the log10(cutoff_ratio) as the independent variable.
//! Tabulated values are stored and interpolated as 64-bit floats, and only converted
//! to 32-bit floats at the final stage of calculation.
//!
//! Filter coefficients are extracted from scipy's state-space representations,
//! which are the result of a bilinear transform of the transfer function polynomials.
//!
//! | Filter | Min. Cutoff Ratio | Max. Cutoff Ratio |
//! |--------|-------------------|-------------------|
//! | Butter1| 10^-5             | 0.4               |
//! | Butter2| 10^-3             | 0.4               |
//! | Butter3| 10^-2.25 (~0.006) | 0.4               |
//! | Butter4| 10^-1.5 (~0.032)  | 0.4               |
//! | Butter5| 10^-1.5 (~0.032)  | 0.4               |
//! | Butter6| 10^-1.25 (~0.06)  | 0.4               |
//!
//! # License
//! Licensed under either of
//!
//! - Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
//! - MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
//!
//! at your option.

#![cfg_attr(not(feature = "std"), no_std)]
mod generated;
pub use generated::butter::butter1::butter1;
pub use generated::butter::butter2::butter2;
pub use generated::butter::butter3::butter3;
pub use generated::butter::butter4::butter4;
pub use generated::butter::butter5::butter5;
pub use generated::butter::butter6::butter6;

// Required for float conversions from 64 to 32 bit and for f32::log10
#[cfg(not(feature = "std"))]
use num_traits::Float;

/// `std` is required for tests, but is not a default feature.
/// To allow the library to compile with default features,
/// tests that require `std` are feature-gated.
/// This test makes sure we do not skip the real tests.
#[cfg(test)]
#[cfg(not(feature = "std"))]
mod test {
    #[test]
    fn require_std_for_tests() {
        panic!("`std` feature is required for tests")
    }
}

/// A simple array with large memory alignment because it will be accessed
/// often in a loop
#[repr(align(32))]
struct AlignedArray<const N: usize>([f32; N]);

/// Single-Input-Single-Output, Infinite Impulse Response filter,
/// normalized to a sample time interval of 1.0
#[repr(align(32))]
pub struct SisoIirFilter<const ORDER: usize> {
    // Aligning the struct will usually keep these scalar fields aligned well enough
    /// Latest state estimate
    y: f32,
    /// State-space `D` scalar
    d: f32,

    /// Internal state buffer, state-space `X(k)` or `X(k-1)`
    x: AlignedArray<ORDER>,
    /// Nontrivial row of state-space `A` matrix in canonical form
    a: AlignedArray<ORDER>,
    /// State-space `C` vector
    c: AlignedArray<ORDER>,
}

impl<const ORDER: usize> SisoIirFilter<ORDER> {
    /// Evaluate the next estimated value based on the latest measurement
    /// in 4N+1 floating-point ops for a filter of order N.
    #[inline]
    pub fn update(&mut self, u: f32) -> f32 {
        // Both multiply-and-sum loops could be turned into chained mul-add,
        // but microcontrollers mostly don't have
        // mul-add instructions as of the year 2024, so using mul_add here
        // would cause severe performance regression. Using chained mul-add
        // on an unknown number of points also requires a recursion that
        // may not be tail-call optimized because it contains a branch,
        // resulting in further increased stack usage and reduced throughput.

        // Y(k) = CX(k-1) + DU(k)
        // 2N+1 float ops
        // Sum starting with d*u because this term is the smallest,
        // and `c` terms are ordered from smallest to largest.
        // Summation from smallest to largest terms improves
        // float roundoff error.
        self.y = self
            .c
            .0
            .iter()
            .zip(self.x.0.iter())
            .fold(self.d * u, |acc, (cval, xval)| acc + cval * xval);

        // X(k) = AX(k-1) + BU(k)
        // `B` in canonical form is like [1, 0, ...] and just selects the one nonzero value in `U`,
        // which is the latest raw measurement.
        // 2N float ops
        let x0 = self
            .a
            .0
            .iter()
            .rev()
            .zip(self.x.0.iter().rev())
            .map(|(aval, xval)| aval * xval)
            .sum::<f32>()
            + u; // The one nontrivial element from BU(k)
        (0..ORDER - 1).for_each(|i| self.x.0[ORDER - 1 - i] = self.x.0[ORDER - 2 - i]); // Time delays
        self.x.0[0] = x0;

        self.y
    }

    pub fn new(a: &[f32], c: &[f32], d: f32) -> Self {
        let mut a_ = [0.0; ORDER];
        a_.copy_from_slice(a);

        let mut c_ = [0.0; ORDER];
        c_.copy_from_slice(c);

        Self {
            y: 0.0,
            x: AlignedArray([0.0; ORDER]),
            a: AlignedArray(a_),
            c: AlignedArray(c_),
            d,
        }
    }

    /// Build a new low-pass with coefficients interpolated on baked tables.
    /// After interpolation, the `C` vector is scaled by a (hopefully) small
    /// amount to more closely produce unity steady-state gain.
    pub fn new_interpolated(
        cutoff_ratio: f64,
        log10_cutoff_ratio_grid: &[f64],
        avals: &[&[f64]],
        cvals: &[&[f64]],
        dvals: &[f64],
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

        let mut a = [0.0; ORDER];
        let mut c = [0.0; ORDER];

        // Interpolate `A` values
        for i in 0..ORDER {
            a[i] = interpn::MulticubicRectilinear::<'_, _, 1>::new(
                &[log10_cutoff_ratio_grid],
                avals[i],
                true,
            )?
            .interp_one(&[log10_cutoff_ratio])? as f32;
        }

        // Interpolate `C` values
        for i in 0..ORDER {
            c[i] = interpn::MulticubicRectilinear::<'_, _, 1>::new(
                &[log10_cutoff_ratio_grid],
                cvals[i],
                true,
            )?
            .interp_one(&[log10_cutoff_ratio])? as f32;
        }

        // Interpolate `D` value
        let d = interpn::MulticubicRectilinear::<'_, _, 1>::new(
            &[log10_cutoff_ratio_grid],
            dvals,
            true,
        )?
        .interp_one(&[log10_cutoff_ratio])? as f32;

        // Scale `C` to enforce unity gain at zero frequency, that is,
        // that the step response should converge exactly.
        //
        // First, find scalar `xs` for every entry in filter state vector `x`
        // such that x(k) == x(k-1).
        let asum = a.iter().sum::<f32>() as f64;
        let xs = 1.0 / (1.0 - asum);
        // Then, find what the sum of `C` _should_ be to produce unity gain,
        // and use that value to calculate a scale factor for the existing `C`.
        let csum_desired = (1.0 - d as f64) / xs;
        let csum = c.iter().sum::<f32>() as f64;
        let scale_factor = (csum_desired / csum) as f32;
        // Finally, scale `C` to, as closely as possible, produce unity gain.
        c.iter_mut().for_each(|v| *v *= scale_factor);

        Ok(Self::new(&a, &c, d))
    }
}
