//! Infinite-impulse-response filtering on 32-bit floats

// Required for float conversions from 64 to 32 bit and for f32::log10
#[cfg(not(feature = "std"))]
#[allow(unused_imports)]
use num_traits::Float;

use super::{AlignedArray, Ring};

/// Single-Input-Single-Output, Infinite Impulse Response filter,
/// normalized to a sample time interval of 1.0
#[derive(Clone, Copy, Default)]
pub struct SisoIirFilter<const ORDER: usize> {
    // Aligning the struct will usually keep these scalar fields aligned well enough
    /// Latest state estimate
    y: f32,
    /// State-space `D` scalar
    d: f32,

    /// Internal state buffer, state-space `X` vector
    x: Ring<f32, ORDER>,
    /// Nontrivial row of state-space `A` matrix in canonical form
    a: AlignedArray<f32, ORDER>,
    /// State-space `C` vector
    c: AlignedArray<f32, ORDER>,
}

impl<const ORDER: usize> SisoIirFilter<ORDER> {
    /// Evaluate the next estimated value based on the latest measurement
    /// in 4N-1 floating-point ops for a filter of order N.
    #[inline]
    pub fn update(&mut self, u: f32) -> f32 {
        // Y(k) = CX(k-1) + DU(k)
        // 2N float ops
        // Sum starting with d*u because this term is the smallest,
        // and `c` terms are ordered from smallest to largest.
        // Summation from smallest to largest terms improves
        // float roundoff error.
        self.y = self.c.dot(&self.x, self.d * u);

        // X(k) = AX(k-1) + BU(k)
        // `B` in canonical form is like [1, 0, ...] and just selects the one nonzero value in `U`,
        // which is the latest raw measurement.
        // 2N-1 float ops
        let x0 = self.a.dot(&self.x, u);

        // Update x0 and apply time delays
        self.x.push(x0);

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
            x: Ring::new(0.0),
            a: AlignedArray(a_),
            c: AlignedArray(c_),
            d,
        }
    }

    /// Build a new low-pass with coefficients interpolated on baked tables.
    /// After interpolation, the `C` vector is scaled by a (hopefully) small
    /// amount to more closely produce unity steady-state gain.
    ///
    /// Note that this procedure will not work for high-pass or band-pass filters,
    /// which are incompatible with the error-correction procedure.
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
            .interp_one([log10_cutoff_ratio])? as f32;
        }

        // Interpolate `C` values
        for i in 0..ORDER {
            c[i] = interpn::MulticubicRectilinear::<'_, _, 1>::new(
                &[log10_cutoff_ratio_grid],
                cvals[i],
                true,
            )?
            .interp_one([log10_cutoff_ratio])? as f32;
        }

        // Interpolate `D` value
        let d = interpn::MulticubicRectilinear::<'_, _, 1>::new(
            &[log10_cutoff_ratio_grid],
            dvals,
            true,
        )?
        .interp_one([log10_cutoff_ratio])? as f32;

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

        // Reverse order to speed up later evaluations
        a.reverse();
        c.reverse();

        Ok(Self::new(&a, &c, d))
    }

    /// Set filter internal state to the steady value
    /// achieved for input `u`. For filters with unity steady-state gain,
    /// this will also produce an output reading of `u`.
    pub fn set_steady_state(&mut self, u: f32) {
        let asum: f32 = self.a.0.iter().sum();
        let xs = u / (1.0 - asum); // Scalar constant value for `X` entries
        self.x = Ring::new(xs);
    }

    /// Latest output
    pub fn y(&self) -> f32 {
        self.y
    }

    /// Read-only access to internal state
    pub fn x(&self) -> &Ring<f32, ORDER> {
        &self.x
    }

    /// Read-only access to nontrivial row of state-space `A`
    pub fn a(&self) -> &AlignedArray<f32, ORDER> {
        &self.a
    }

    /// Read-only access to nontrivial row of state-space `C`
    pub fn c(&self) -> &AlignedArray<f32, ORDER> {
        &self.c
    }

    /// Read-only access to nontrivial entry in state-space `D`
    pub fn d(&self) -> f32 {
        self.d
    }
}
