#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
pub mod fir;
pub mod fractional_delay;
pub mod iir;
pub mod median;

use crunchy::unroll;
pub use fir::SisoFirFilter;
pub use fractional_delay::polynomial_fractional_delay;
pub use iir::SisoIirFilter;
pub use median::MedianFilter;

pub mod generated;
pub use generated::butter::butter1::butter1;
pub use generated::butter::butter2::butter2;
pub use generated::butter::butter3::butter3;
pub use generated::butter::butter4::butter4;
pub use generated::butter::butter5::butter5;
pub use generated::butter::butter6::butter6;

use num_traits::Num;

/// A simple array with large memory alignment because it will be accessed
/// often in a loop, with methods specialized for filter evaluation.
#[derive(Clone, Copy)]
#[repr(align(8))]
pub struct AlignedArray<T, const N: usize>([T; N]);

impl<T: Copy + Num, const N: usize> AlignedArray<T, N> {
    /// Multiply-and-sum between this array and a target ring buffer, starting
    /// with the most recent sample and the first element of this array
    /// and finishing with the least recent sample and the last element of this array.
    ///
    /// A starting value can be provided for the summation,
    /// which can be helpful for fine-tuning floating-point error.
    #[inline]
    pub fn dot(&self, buf: &Ring<T, N>, start: T) -> T {
        const {assert!(N < 128, "N > 128 not supported")}

        let other = buf.buf();
        let mut acc = start;
        unroll! {
            for i < 128 in 0..N {
                // Iterate in reverse for improved numerics (sum smallest values first)
                let j = (N-1) - i;
                acc = acc + self.0[j] * other[j];
            }
        }

        acc
    }
}

/// Ring buffer.
/// Most recent sample is stored first.
#[derive(Clone, Copy)]
#[repr(transparent)]
pub struct Ring<T: Copy, const N: usize> {
    buf: AlignedArray<T, N>,
}

impl<T: Copy, const N: usize> Ring<T, N> {
    /// Initialize with buffer populated with constant value
    pub fn new(value: T) -> Self {
        Self {
            buf: AlignedArray([value; N])
        }
    }

    /// Replace the oldest value in the buffer with a new value
    pub fn push(&mut self, value: T) {
        unroll! {
            for i < 128 in 1..N {
                self.buf.0[N-i] = self.buf.0[N-i-1];
            }
        }
        self.buf.0[0] = value;
    }

    /// The whole internal buffer, with no indication of the current index
    fn buf(&self) -> &[T; N] {
        &self.buf.0
    }
}

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
