# Changelog

## 0.6.0 - 2025-11-29

## Added
* `sos::SisoSosFilter`, an implementation of cascaded second order sections
  * includes `new_interpolated` initializer, which interpolates the SOS coefficients versus cutoff ratio from a lookup table
  * generic over `f32` and `f64`
  * FMA and non-FMA implementations, configured by the existing `fma` feature
* functions to generate SOS Butterworth lowpass filters from lookup tables: `sos::{butter2, butter4, butter6}`
  * provides orders 2, 4 and 6 for `f32` and `f64`
  * `cutoff_ratio` region of validity is enforced and tested for each combination of filter order and float type
  * `cutoff_ratio` region of validity is larger than existing Butterworth implementation for orders 4 and up
    * e.g. for order 4, `f32`: old 0.032, new 0.005
  * SOS *does not* improve region of validity for order 2 filter
* python script to autogenerate the SOS Butterworth lookup tables: `scripts/generate_butter_sos_tables.py`
* benchmarking with criterion, and a benchmark of `sos::butter4::<f64>`

## Changed

* !Renamed `flaw::SisoIirFilter::initialize` to `set_steady_state`

## 0.5.0 - 2025-10-19

Using FMA to accelerate and reduce roundoff error in dot products.

### Added

* Add `fma` feature that enables the use of fused multiply-add for dot products
  * Disabled by default because this causes a severe performance regression unless FMA is enabled for the target

### Changed

* !Add trait bound on `MulAdd<Output = T>` for numeric types in FIR and IIR filters
* !Make `std` a default feature to enable publishing and docs to work properly

## 0.4.0 - 2025-10-03

### Changed

* !Add bound on Num trait for generic fields
* Implement Default on all structs
* Update deps

Various performance-related improvements.

## 0.3.0 - 2025-08-31

Various performance-related improvements.
~20% improvement in speed based on testing on an STM32H7 microcontroller.

### Changed

* !Store filter coeffs and taps in opposite order (most-recent-last) to improve vectorization
* !Const-unroll compute intensive loops
    * !Limit size of FIR filters to <128 taps
    * !Limit order of polynomial fractional delay filters to <=10
* Add `u` value as initial for A*x dot product
    * (Potentially substantially) Improves float error for small values of latest measurement
      `u` at the expense of slightly worse float error for large `u`
* Update flop per eval count to 4N-1
* Update rust edition to 2024

## 0.2.5 - 2025-05-25

### Changed

* Make `Ring` and `AlignedArray` pub to support read-only access to filter internals

### Added

* Add FIR module with generic number field type
* Add Lagrange interpolator fractional delay filter construction
* Add functions for read-only access to IIR filter internals

## 0.2.4 - 2025-02-06

### Changed

* Change large const arrays to static to avoid inlining excessively large data

## 0.2.3 - 2025-01-01

### Changed

* Use ring buffer for IIR and median filter state storage
* Lower align of SisoIirFilter struct, since the internal storage is already aligned

### Added

* Implement N-point median filter for odd N >= 3

## 0.2.2 - 2024-12-27

### Changed

* Derive Clone and Copy on `SisoIirFilter` to allow faster construction of batches of identical filters
* Make `generated` module pub to provide access to min/max cutoff ratio bounds for generated filters
* Loosen version requirement on num_traits dep

## 0.2.1 - 2024-12-26

### Changed

* Lower align of `SisoIirFilter` struct and internal state arrays to 8

### Added

* Add `SisoIirFilter.initialize(u)` to set the internal state vector to the steady-state values associated with a given measurement `u`

## 0.2.0 - 2024-07-10

### Changed

* Implement rescaling of `C` vector during initialization of interpolated filters to correct gain for interpolation error
* Tighten tolerance on gain of generated butterworth filters to 0.01% at minimum cutoff ratio
* Decrease range of validity for generated butterworth filters

## 0.1.0 - 2024-04-07

### Added

* Initial release
