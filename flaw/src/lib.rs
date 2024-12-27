#![doc = include_str!("../README.md")]
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
#[repr(align(8))]
struct AlignedArray<const N: usize>([f32; N]);

/// Single-Input-Single-Output, Infinite Impulse Response filter,
/// normalized to a sample time interval of 1.0
#[repr(align(8))]
pub struct SisoIirFilter<const ORDER: usize> {
    // Aligning the struct will usually keep these scalar fields aligned well enough
    /// Latest state estimate
    y: f32,
    /// State-space `D` scalar
    d: f32,

    /// Internal state buffer, state-space `X` vector
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

    /// Populate a new filter with arbitrary state-space `A`, `C`, and `D`
    /// in canonical form s.t. `A` values are the top row only
    /// and `B` is assumed to be unity.
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

    /// Initialize filter internal state to the steady value
    /// achieved for input `u`. For filters with unity steady-state gain,
    /// this will also produce an output reading of `u`.
    pub fn initialize(&mut self, u: f32) {
        let asum: f32 = self.a.0.iter().sum();
        let xs = u / (1.0 - asum); // Scalar constant value for `X` entries
        self.x.0.iter_mut().for_each(|x| *x = xs);
    }
}

#[cfg(feature = "std")]
#[cfg(test)]
mod test {
    /// Test initialization to a given input value
    #[test]
    fn test_initialize() {
        let e = core::f64::consts::E as f32;
        let mut f = super::butter2(0.2).unwrap();
        f.initialize(e);
        assert!((e - f.update(e)).abs() / e < 1e-6);
    }
}
