# Chanelog

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
