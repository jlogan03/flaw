#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
pub mod fir;
pub mod fractional_delay;
pub mod iir;
pub mod median;
pub mod sos;

use core::ops::{Index, IndexMut};
use crunchy::unroll;
pub use fir::SisoFirFilter;
pub use fractional_delay::polynomial_fractional_delay;
pub use iir::SisoIirFilter;
pub use median::MedianFilter;
pub use sos::SisoSosFilter;

pub mod generated;
pub use generated::butter::butter1::butter1;
pub use generated::butter::butter2::butter2;
pub use generated::butter::butter3::butter3;
pub use generated::butter::butter4::butter4;
pub use generated::butter::butter5::butter5;
pub use generated::butter::butter6::butter6;

use num_traits::{MulAdd, Num};

/// A simple array with large memory alignment because it will be accessed
/// often in a loop, with methods specialized for filter evaluation.
#[derive(Clone, Copy, Debug)]
#[repr(align(8))]
pub struct AlignedArray<T, const N: usize>([T; N]);

impl<T: Copy + Num, const N: usize> Default for AlignedArray<T, N> {
    fn default() -> Self {
        Self([T::zero(); N])
    }
}

impl<T: Copy + Num + MulAdd<Output = T>, const N: usize> AlignedArray<T, N> {
    /// Multiply-and-sum between this array and a target ring buffer.
    ///
    /// A starting value can be provided for the summation,
    /// which can be helpful for fine-tuning floating-point error.
    #[inline]
    pub fn dot(&self, buf: &Ring<T, N>, start: T) -> T {
        const { assert!(N < 129, "N > 128 not supported") }

        // Multiply-and-sum loops could be turned into chained mul-add,
        // but microcontrollers mostly don't have mul-add instructions
        // as of the year 2025, so using mul_add here would cause severe
        // performance regression.
        let other = buf.buf();
        let mut acc = start;

        #[cfg(not(feature = "fma"))]
        unroll! {
            for i < 128 in 0..N {
                acc = acc + self.0[i] * other[i];
            }
        }

        #[cfg(feature = "fma")]
        unroll! {
            for i < 128 in 0..N {
                acc = self.0[i].mul_add(other[i], acc);
            }
        }

        acc
    }
}

impl<T, const N: usize> Index<usize> for AlignedArray<T, N> {
    type Output = T;
    #[inline]
    fn index(&self, idx: usize) -> &Self::Output {
        &self.0[idx]
    }
}

impl<T, const N: usize> IndexMut<usize> for AlignedArray<T, N> {
    #[inline]
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        &mut self.0[idx]
    }
}

/// Ring buffer.
/// Most recent sample is stored last.
#[derive(Clone, Copy, Default, Debug)]
#[repr(transparent)]
pub struct Ring<T: Num + Copy, const N: usize> {
    buf: AlignedArray<T, N>,
}

impl<T: Num + Copy, const N: usize> Ring<T, N> {
    /// Initialize with buffer populated with constant value
    pub fn new(value: T) -> Self {
        Self {
            buf: AlignedArray([value; N]),
        }
    }

    /// Replace the oldest value in the buffer with a new value.
    /// Despite the unnecessary copies compared to a true ring buffer, this is
    /// faster overall in filtering applications where one dot product is performed
    /// per push, and the cost of splitting the dot product into two segments is large
    /// compared to the savings from avoiding copies.
    pub fn push(&mut self, value: T) {
        unroll! {
            for i < 128 in 0..(N-1) {
                self.buf.0[i] = self.buf.0[i + 1];
            }
        }
        self.buf.0[N - 1] = value;
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
        f.set_steady_state(e);
        assert!((e - f.update(e)).abs() / e < 1e-6);
    }
}
