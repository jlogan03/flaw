#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
pub mod iir;
pub mod fir;
pub mod median;
pub mod fractional_delay;

pub use iir::SisoIirFilter;
pub use fir::SisoFirFilter;
pub use median::MedianFilter;
pub use fractional_delay::polynomial_fractional_delay;

pub mod generated;
pub use generated::butter::butter1::butter1;
pub use generated::butter::butter2::butter2;
pub use generated::butter::butter3::butter3;
pub use generated::butter::butter4::butter4;
pub use generated::butter::butter5::butter5;
pub use generated::butter::butter6::butter6;

use num_traits::Num;

/// A simple array with large memory alignment because it will be accessed
/// often in a loop
#[derive(Clone, Copy)]
#[repr(align(8))]
struct AlignedArray<T, const N: usize>([T; N]);

impl<T: Copy + Num, const N: usize> AlignedArray<T, N> {
    /// Multiply-and-sum between this array and a target ring buffer, starting
    /// with the most recent sample and the first element of this array
    /// and finishing with the least recent sample and the last element of this array.
    /// 
    /// A starting value can be provided, which can be helpful for fine-tuning floating-point error.
    #[inline]
    pub fn dot(&self, buf: &Ring<T, N>, start: T) -> T {
        // Split buffers into compatible slices
        let x_parts = buf.buf_parts(); // [first, second] both reversed
        let n = x_parts.0.len(); // ring buffer split point w.r.t. contiguous vectors
        let a_parts = self.0.split_at(n);

        let mut out = start;
        //    Sum the second contiguous segment first because it will have smaller values
        //    for low-pass IIR filters, so this substantially improves float precision.
        //    The tests don't pass without this as of 2025-05-14!
        a_parts
            .1
            .iter()
            .zip(x_parts.1.iter().rev())
            .for_each(|(&aval, &xval)| out = out + aval * xval);
        //    Sum the first contiguous segment  
        a_parts
            .0
            .iter()
            .zip(x_parts.0.iter().rev())
            .for_each(|(&aval, &xval)| out = out + aval * xval);

        out
    }
}

/// Ring buffer
#[derive(Clone, Copy)]
struct Ring<T: Copy, const N: usize> {
    buf: AlignedArray<T, N>,
    next: usize,
}

impl<T: Copy, const N: usize> Ring<T, N> {
    /// Initialize with buffer populated with constant value
    fn new(value: T) -> Self {
        Self {
            buf: AlignedArray([value; N]),
            next: 0,
        }
    }

    /// Replace the oldest value in the buffer with a new value
    fn push(&mut self, value: T) {
        self.buf.0[self.next] = value;
        self.next += 1;
        if self.next == N {
            self.next = 0;
        }
    }

    /// The whole internal buffer, with no indication of the current index
    fn buf(&self) -> &[T; N] {
        &self.buf.0
    }

    /// Contiguous parts of the buffer split across the next insertion point
    /// s.t. the first slice includes values from index 0 to the most recently
    /// inserted value.
    ///
    /// To iterate over the items in order of insertion from most recent to oldest,
    /// loop over the first slice, then the second, both in reverse.
    fn buf_parts(&self) -> (&[T], &[T]) {
        self.buf.0.split_at(self.next)
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
