"""Generation of Butterworth filter coefficient tables"""
import os
from pathlib import Path

from scipy.signal import butter, dlti, StateSpace, dlsim
import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).parent

n = 100
orders = [1, 2, 3, 4, 5, 6]
min_log10_cutoffs = [-5, -3, -2.25, -1.5, -1.5, -1.25]
max_log10_cutoff = float(np.log10(0.4))

f_ref = 10.0**max_log10_cutoff  # [dimensionless] reference frequency ratio for cutoff testing
tmax = 100.0  # Brute force method for getting good resolution on attenuation
cutoff_test_input = np.sin(2.0 * np.pi * f_ref * np.linspace(0.0, tmax, int(tmax) + 1, endpoint=True))

modfpath = here / f"../flaw/src/generated/butter/mod.rs"
with open(modfpath, "w") as f:
    f.write("//! Butterworth filters based on baked filter coefficient tables.\n//! This file is autogenerated.\n\n")
    for order in orders:
        f.write(f"pub mod butter{order};\n")

for order, min_log10_cutoff in zip(orders, min_log10_cutoffs):
    # Cutoff frequency as a fraction of sampling frequency.
    # All values must be less than (not less-than-or-equal) the Nyquist frequency (half of sampling),
    # hence the final division.
    # Higher-order filters see some numerical difficulty in producing the coeffs, and have their domain trimmed
    # to avoid including junk outputs.
    cutoff_ratios = np.logspace(min_log10_cutoff, max_log10_cutoff, n, endpoint=True)  # [dimensionless]

    # Generate discrete LTI systems with a samplerate of 1
    # (in any units, because the coeffs are generated relative to that frequency).
    dltis = [dlti(*butter(N=order, Wn=c, fs=1.0), dt=1.0) for c in cutoff_ratios]
    systems = [StateSpace(s) for s in dltis]
    amats = [s.A for s in systems]
    bmats = [s.B for s in systems]
    cmats = [s.C for s in systems]
    dmats = [s.D for s in systems]

    # Generate test output against a subset of the tabulated cutoff ratios
    _t, cutoff_test_output = dlsim(dlti(*butter(N=order, Wn=f_ref, fs=1.0), dt=1.0), cutoff_test_input)
    nstepmin = int(1/cutoff_ratios[0]) * 10  # Make sure it is converged
    nstepmax = int(1/cutoff_ratios[-1])
    _t, step_test_min = dlsim(dlti(*butter(N=order, Wn=cutoff_ratios[0], fs=1.0), dt=1.0), np.ones(nstepmin))
    _t, step_test_max = dlsim(dlti(*butter(N=order, Wn=cutoff_ratios[-1], fs=1.0), dt=1.0), np.ones(nstepmax))

    #
    # Bake static tables
    #
    outfile = here / f"../flaw/src/generated/butter/butter{order}.rs"

    with open(outfile, "w") as f:
        f.write(f"//! Butterworth filter of order {order}.\n")
        f.write(f"//! Region of validity: cutoff ratio from {cutoff_ratios[0]:.2e} to {cutoff_ratios[-1]:.2e} .\n")
        f.write(f"//! This file is autogenerated.\n")
        f.write("#![allow(clippy::style)]\n\n")
        f.write("use crate::SisoIirFilter;\n\n")

        f.write(f"/// Minimum tabulated cutoff ratio\n")
        f.write("#[allow(dead_code)]\n")
        f.write(f"pub const MIN_CUTOFF_RATIO: f64 = {float(cutoff_ratios[0])};\n\n")

        f.write(f"/// Maximum tabulated cutoff ratio\n")
        f.write("#[allow(dead_code)]\n")
        f.write(f"pub const MAX_CUTOFF_RATIO: f64 = {float(cutoff_ratios[-1])};\n\n")

        f.write(f"/// Initialise a Butterworth filter of order {order} by interpolating the coefficients from stored tables.\n")
        f.write(f"/// Cutoff ratio is the dimensionless ratio of the cutoff frequency to the sampling frequency.\n")
        f.write(f"/// Region of validity: cutoff ratio from {cutoff_ratios[0]:.2e} to {cutoff_ratios[-1]:.2e}\n")
        f.write(f"pub fn butter{order}(cutoff_ratio: f64) -> Result<SisoIirFilter<{order}>, &'static str>" " {\n")
        avals_string = "".join(["    let avals = &["] + [f"&AVALS[{i}][..], " for i in range(order)] + ["];\n"])
        cvals_string = "".join(["    let cvals = &["] + [f"&CVALS[{i}][..], " for i in range(order)] + ["];\n"])
        f.write(avals_string)
        f.write(cvals_string)
        f.write("    SisoIirFilter::new_interpolated(cutoff_ratio, &LOG10_CUTOFF_RATIOS, avals, cvals, &DVALS)\n")
        f.write("}\n\n")

        logcrlist = [np.log10(x) for x in cutoff_ratios]
        f.write("/// [dimensionless] Log base-10 of cutoff ratios, to improve float precision during interpolation\n")
        f.write("#[rustfmt::skip]\n")
        f.write(f"const LOG10_CUTOFF_RATIOS: [f64; {n}] = {logcrlist};\n\n")

        dvals = [d[0][0] for d in dmats]
        f.write("/// State-Space `D` 1x1 matrix\n")
        f.write("#[rustfmt::skip]\n")
        f.write(f"const DVALS: [f64; {n}] = {dvals};\n\n")

        avals = [[a[0,i] for a in amats] for i in range(order)]
        f.write(f"/// State-Space `A` matrix, first row\n")
        f.write("#[rustfmt::skip]\n")
        f.write(f"const AVALS: [[f64; {n}]; {order}] = {avals};\n\n")

        cvals = [[c[0,i] for c in cmats] for i in range(order)]
        f.write(f"/// State-Space `C` vector\n")
        f.write("#[rustfmt::skip]\n")
        f.write(f"const CVALS: [[f64; {n}]; {order}] = {cvals};\n\n")

        f.write('#[cfg(feature = "std")]\n')
        f.write("#[cfg(test)]\n")
        f.write("#[rustfmt::skip]\n")
        f.write("mod test {\n")
        f.write("    use super::*;\n")
        f.write(f"    const CUTOFF_TEST_INPUT: [f32; {cutoff_test_input.size}] = {[x for x in cutoff_test_input]};\n")
        f.write(f"    const CUTOFF_TEST_OUTPUT: [f32; {cutoff_test_output.size}] = {[x for x in cutoff_test_output.T[0]]};\n")
        f.write(f"    const STEP_TEST_MIN_OUTPUT: f32 = {float(step_test_min[-1][0])};\n")
        f.write(f"    const STEP_TEST_MAX_OUTPUT: f32 = {float(step_test_max[-1][0])};\n\n")
        f.write("    #[test]\n")
        f.write("    fn test() {\n")
        f.write(f'        println!("order {order}");\n')
        f.write(f"        let mut filter = butter{order}({f_ref}).unwrap();\n")
        f.write("        let out = (0..CUTOFF_TEST_INPUT.len()).map(|i| {filter.update(CUTOFF_TEST_INPUT[i])}).collect::<Vec<f32>>();\n")
        f.write("        // Check overall match to reference output to catch phase error, etc\n")
        f.write("        (0..CUTOFF_TEST_INPUT.len()).for_each(|i| { let expected = CUTOFF_TEST_OUTPUT[i]; let rel_err = (out[i] - expected).abs() / expected.abs().max(1e-4); assert!(rel_err < 0.05); });\n")
        f.write("        // Check approximate attenuation at cutoff frequency; should be -3dB or 1/sqrt(2) magnitude\n")
        f.write("        let maxmag = out.iter().fold(0.0_f32, |a, b| a.abs().max(b.abs()));\n")
        f.write("        let attenuation_rel_err = (maxmag - 0.707).abs() / 0.707;\n")
        f.write('        println!("attenuation rel err {attenuation_rel_err}");\n')
        f.write("        assert!(attenuation_rel_err < 0.05);\n")
        f.write("        // Check convergence of step responses at min and max tabulated cutoff\n")
        f.write(f"        let mut filtermin = butter{order}(MIN_CUTOFF_RATIO).unwrap();\n")
        f.write(f"        (0..{nstepmin - 1})"".for_each(|_| {filtermin.update(1.0);});\n")
        f.write("        let step_min_final = filtermin.update(1.0);\n")
        f.write("        let step_min_rel_err = (step_min_final - STEP_TEST_MIN_OUTPUT).abs() / STEP_TEST_MIN_OUTPUT;\n")
        f.write('        println!("step min rel err {step_min_rel_err}");\n')
        f.write("        assert!(step_min_rel_err < 1e-3);\n")
        f.write(f"        let mut filtermax = butter{order}(MAX_CUTOFF_RATIO).unwrap();\n")
        f.write(f"        (0..{nstepmax - 1})"".for_each(|_| {filtermax.update(1.0);});\n")
        f.write("        let step_max_final = filtermax.update(1.0);\n")
        f.write("        let step_max_rel_err = (step_max_final - STEP_TEST_MAX_OUTPUT).abs() / STEP_TEST_MAX_OUTPUT;\n")
        f.write('        println!("step max rel err {step_max_rel_err}");\n')
        f.write("        assert!(step_max_rel_err < 1e-6);\n")
        f.write("    }\n")
        f.write("}\n")

    #
    # Plot
    #
    fig, axes = plt.subplots(1, 3, sharex=True, figsize=(12,5))
    # Plot top row of A, which is in canonical form s.t. this
    # is the only nontrivial row and the others are all time delays
    plt.sca(axes[0])
    for i in range(order):
        plt.semilogx(cutoff_ratios, [a[0,i] for a in amats])
    plt.title("A")
    plt.xlabel("Normalized Cutoff Frequency W/Ws [dimensionless]")

    # C is a fully populated vector
    plt.sca(axes[1])
    for i in range(order):
        plt.loglog(cutoff_ratios, [c[0, i] for c in cmats])
    plt.title("C")
    plt.xlabel("Normalized Cutoff Frequency W/Ws [dimensionless]")

    # D is a scalar in a 2D matrix type
    plt.sca(axes[2])
    plt.loglog(cutoff_ratios, [d[0][0] for d in dmats])
    plt.title("D")
    plt.xlabel("Normalized Cutoff Frequency W/Ws [dimensionless]")

    plt.suptitle(f"Discrete Butterworth, Order = {order}")
    plt.tight_layout(rect=[0, 0.03, 1, 0.99])

plt.show()