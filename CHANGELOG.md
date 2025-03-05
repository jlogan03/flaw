# Changelog

## 0.2.5 - 2025-03-05

## Changed

* Update gain error correction procedure to reduce float error
* Slightly widen butter2 minimum cutoff
* Reduce number of table entries to 30 from 100

## Added

* Add `StagedSisoIirFilter` struct to initialize and update multiple identical filter stages
* Add `butter{N}_2stage` variants that create a staged filter with the desired aggregate cutoff and combined order of `2*{N}`
* Add `butter{N}::LOG10_ROOT_CUTOFF_RATIO` statics describing the normalized frequency where the gain is sqrt(cutoff) for a given corresponding cutoff in `butter{N}::LOG10_CUTOFF_RATIO`, which allows looking up the appropriate cutoff ratio for individual filter stages when building multi-stage filters

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
