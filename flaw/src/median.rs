//! Median filter on arbitrary partially-ordered type.

use super::Ring;

/// In-place insertion sort method.
///
/// Yoinked without modification from the `sorts` crate
/// in order to use in a no-std environment which is not
/// supported by the original.
fn insertion_sort<T: PartialOrd>(s: &mut [T]) {
    for i in 1..s.len() {
        let mut j = i;
        while j > 0 && s[j - 1] > s[j] {
            s.swap(j - 1, j);
            j -= 1;
        }
    }
}

/// Simple median filter on arbitrary ordered type.
/// Requires an odd number of points to make the median unique.
///
/// Due to the fact that rust's slice::sort method is part of `alloc`,
/// this uses an insertion sort method instead. Because of this, it is
/// not recommended for large `N`, but is likely faster than other sort methods
/// when used for appropriately small `N`.
///
/// Note that while this method only requires `PartialOrd`, which allows it
/// to be used with float values, care must be taken not to poison the buffer
/// with NaN values, which will cause the identification of the median to fail
/// silently. Because the .is_finite() method is only available for floats,
/// and not for other numeric types, guards against this failure must be
/// implemented on the inputs given to this filter, and can't be implemented
/// inside the filter update.
pub struct MedianFilter<T: PartialOrd + Copy, const N: usize> {
    vals: Ring<T, N>,
    buf: [T; N],
    imid: usize,
}

impl<T: PartialOrd + Copy, const N: usize> MedianFilter<T, N> {
    /// Initialize with the buffer fully populated with the supplied value `v`
    pub fn new(v: T) -> Self {
        const {
            assert!(N >= 3);
            assert!(N % 2 == 1);
        }

        let imid = N / 2; // Index of middle of buffer

        Self {
            vals: Ring::new(v),
            buf: [v; N],
            imid,
        }
    }

    /// Push a new value and return the new median
    pub fn update(&mut self, v: T) -> T {
        // Roll input buffer
        self.vals.push(v);

        // Copy, sort, and take the middle element
        // As long as N is small, this will be fast
        self.buf.copy_from_slice(self.vals.buf());
        insertion_sort(&mut self.buf);

        self.buf[self.imid]
    }
}

#[cfg(feature = "std")]
#[cfg(test)]
mod test {
    /// Test initialization to a given input value
    #[test]
    fn test_median() {
        use super::MedianFilter;

        let inp = [2_u32, 5, 6, 7, 8, 4, 9, 0, 1, 3];

        // 3-point
        let expected = [0, 2, 5, 6, 7, 7, 8, 4, 1, 1];
        let mut f = MedianFilter::<u32, 3>::new(0);
        for i in 0..inp.len() {
            let v = f.update(inp[i]);
            let e = expected[i];
            assert!(v == e, "{i}: {v} != {e}");
        }

        // 5-point
        let expected = [0, 0, 2, 5, 6, 6, 7, 7, 4, 3];
        let mut f = MedianFilter::<u32, 5>::new(0);
        for i in 0..inp.len() {
            let v = f.update(inp[i]);
            let e = expected[i];
            assert!(v == e, "{i}: {v} != {e}");
        }
    }
}
