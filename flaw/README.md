# flaw
Embedded signal filtering, no-std and no-alloc compatible.

This library provides a simple method for initializing and updating single-input,
single-output infinite-impulse-response filters using 32-bit floats, as well as
tabulated filter coefficients for some common filters. Filters evaluate in
4N-1 floating-point operations for a filter of order N.

The name `flaw` is short for filter-law, but also refers to the fact that
digital IIR filtering with small floating-point types is an inherently flawed
approach, in that higher-order and lower-cutoff filters produce very small
coefficients that result in floating-point roundoff error. This library mitigates
that problem by providing filter coefficients for a tested
domain of validity. The result is a limited, but useful, range of operation
where these filters can achieve both accuracy and performance as well
as be formulated and initialized in an embedded environment.

## Capabilities

* IIR (f32-only for now)
  * General IIR filter using state-space canonical form
  * Interpolated low-pass filters w/ gain error correction
  * Baked coefficients for Butterworth filters of order 1-6
* FIR (generic number type)
  * General FIR filter
  * Lagrange polynomial fractional-delay filter construction

## Example: Second-Order Butterworth Filter

```rust
// First, choose a cutoff frequency as a fraction of sampling frequency
let cutoff_ratio = 1e-3;

// Construct a filter, interpolating coefficients to that cutoff ratio.
// Initializes internal state to zero by default.
let mut filter = flaw::butter2(cutoff_ratio).unwrap();  // Errors if extrapolating

// Initialize the internal state of the filter
// to match the steady-state associated with some input value.
let initial_steady_measurement = 1.57;  // Some number
filter.initialize(initial_steady_measurement);

// Update the filter with a new raw measurement
let measurement = 0.3145; // Some number
let estimate = filter.update(measurement);  // Latest state estimate
```

## Coefficient Tables

Tabulated filters are tested to enforce

* <0.01% error in converged step response at the minimum cutoff frequency
* <1ppm error in converged step response at the maximum cutoff frequency
* <5% error to -3dB attenuation of a sine input at the cutoff frequency at the maximum cutoff ratio
  * This error appears to be mainly an issue of discretization in test cases, and could be reduced
    by using a better method for testing (fit a sine curve to the result or do gradient-descent
    on a cubic interpolator)

Each filter with tabulated coefficients has a minimum and maximum cutoff ratio.
The minimum value is determined by floating-point error in convergence of a
step response, while the maximum value is determined by the accuracy of attenuation
at the cutoff frequency as the cutoff ratio approaches the Nyquist frequency.

Coefficients for a given filter are interpolated on these tables using a
cubic Hermite method with the log10(cutoff_ratio) as the independent variable.
Tabulated values are stored and interpolated as 64-bit floats, and only converted
to 32-bit floats at the final stage of calculation.

After interpolation, the state-space measurement coefficient vector (`C`) is scaled
to correct steady-state gain for interpolation error, targeting unity gain.

Filter coefficients are extracted from scipy's state-space representations,
which are the result of a bilinear transform of the transfer function polynomials.

| Filter | Min. Cutoff Ratio | Max. Cutoff Ratio |
|--------|-------------------|-------------------|
| Butter1| 10^-4             | 0.4               |
| Butter2| 10^-3             | 0.4               |
| Butter3| 10^-2             | 0.4               |
| Butter4| 10^-1.5 (~0.032)  | 0.4               |
| Butter5| 10^-1.25 (~0.056) | 0.4               |
| Butter6| 0.1               | 0.4               |

## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.
