use num_traits::Num;

use crate::{AlignedArray, Ring};

/// Single-Input-Single-Output, Finite Impulse Response filter.
/// This is simply a running convolution of the taps with the samples.
#[derive(Clone, Copy)]
pub struct SisoFirFilter<const ORDER: usize, T: Num + Copy> {
    /// Internal sample buffer
    x: Ring<T, ORDER>,
    /// Filter taps ordered from most-recent sample to least-recent sample
    taps: AlignedArray<T, ORDER>,
}

impl<const ORDER: usize, T: Num + Copy> SisoFirFilter<ORDER, T> {
    /// Evaluate the next estimated value based on the latest measurement
    /// in 2N-1 floating-point ops for a filter of order N.
    /// This is a simple convolution - push the latest value, then multiply-and-sum.
    #[inline]
    pub fn update(&mut self, u: T) -> T {
        self.x.push(u);
        self.taps.dot(&self.x, T::zero())
    }

    /// Populate a new filter with arbitrary taps.
    /// Filter taps ordered from most-recent sample to least-recent sample.
    pub fn new(taps: &[T]) -> Self {
        let mut taps_ = [T::zero(); ORDER];
        taps_.copy_from_slice(taps);

        Self {
            x: Ring::new(T::zero()),
            taps: AlignedArray(taps_),
        }
    }

    /// Initialize filter internal state to the steady value
    /// achieved for input `u`. For filters with unity steady-state gain,
    /// this will also produce an output reading of `u`.
    pub fn initialize(&mut self, u: T) {
        self.x = Ring::new(u);
    }
}
