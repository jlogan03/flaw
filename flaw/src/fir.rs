use num_traits::Num;

use crate::{AlignedArray, Ring};

/// Single-Input-Single-Output, Finite Impulse Response filter.
/// This is simply a running convolution of the taps with the samples.
#[derive(Clone, Copy)]
pub struct SisoFirFilter<const ORDER: usize, T: Num + Copy> {
    /// Latest output
    y: T,
    /// Internal sample buffer
    x: Ring<T, ORDER>,
    /// Filter taps ordered most-recent-last
    taps: AlignedArray<T, ORDER>,
}

impl<const ORDER: usize, T: Num + Copy> SisoFirFilter<ORDER, T> {
    /// Evaluate the next estimated value based on the latest measurement
    /// in 2N-1 floating-point ops for a filter of order N.
    /// This is a simple convolution - push the latest value, then multiply-and-sum.
    #[inline]
    pub fn update(&mut self, u: T) -> T {
        self.x.push(u);
        self.y = self.taps.dot(&self.x, T::zero());
        self.y
    }

    /// Populate a new filter with arbitrary taps.
    /// Filter taps are ordered most-recent-last.
    pub fn new(taps: &[T]) -> Self {
        let mut taps_ = [T::zero(); ORDER];
        taps_.copy_from_slice(taps);

        Self {
            y: T::zero(),
            x: Ring::new(T::zero()),
            taps: AlignedArray(taps_),
        }
    }

    /// Initialize filter internal state to the steady value
    /// achieved for input `u`. For filters with unity steady-state gain,
    /// this will also produce an output reading of `u`.
    pub fn initialize(&mut self, u: T) {
        self.x = Ring::new(u);
        self.update(u);
    }

    /// Read-only access to taps.
    /// This is, equivalently, the one nontrivial row of either the state-space `A` or `C`
    /// matrix depending on choice of formulation.
    pub fn taps(&self) -> &[T; ORDER] {
        &self.taps.0
    }

    /// Latest output
    pub fn y(&self) -> T {
        self.y
    }

    /// Read-only access to internal state buffer
    pub fn x(&self) -> &Ring<T, ORDER> {
        &self.x
    }
}
