#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
mod iir;
mod median;

pub use iir::{SisoIirFilter, StagedSisoIirFilter};
pub use median::MedianFilter;

pub mod generated;
pub use generated::butter::butter1::butter1;
pub use generated::butter::butter2::butter2;
pub use generated::butter::butter3::butter3;
pub use generated::butter::butter4::butter4;
pub use generated::butter::butter5::butter5;
pub use generated::butter::butter6::butter6;

/// A simple array with large memory alignment because it will be accessed
/// often in a loop
#[derive(Clone, Copy)]
#[repr(align(8))]
struct AlignedArray<T, const N: usize>([T; N]);

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
