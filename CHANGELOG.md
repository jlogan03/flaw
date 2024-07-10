# Chanelog

## 0.2.0 - 2024-07-10

### Changed

* Implement rescaling of `C` vector during initialization of interpolated filters to correct gain for interpolation error
* Tighten tolerance on gain of generated butterworth filters to 0.01% at minimum cutoff ratio
* Decrease range of validity for generated butterworth filters

## 0.1.0 - 2024-04-07

### Added

* Initial release
