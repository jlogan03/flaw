from scipy.signal import butter
import numpy as np
from matplotlib import pyplot as plt

n = 100
sections = 3
order = 2 * sections

# [dimensionless] Cutoff frequency / sample frequency
cutoff_ratios = np.logspace(-4, float(np.log10(0.4)), n, endpoint=True)

# Compute SOS coefficients for each cutoff ratio
sos = np.zeros((sections, 6, n))
for i, c in enumerate(cutoff_ratios):
    sos[:, :, i] = butter(N=order, Wn=c, fs=1.0, btype="low", output="sos")

# TODO generate rust code for these coefficients

# Plot the SOS coefficients
fig, axes = plt.subplots(2, 3, figsize=(10, 8), sharex=True)
for sec in range(sections):
    axes[0, 0].plot(cutoff_ratios, sos[sec, 0, :], label=f"Section {sec + 1}")
    axes[0, 1].plot(cutoff_ratios, sos[sec, 1, :], label=f"Section {sec + 1}")
    axes[0, 2].plot(cutoff_ratios, sos[sec, 2, :], label=f"Section {sec + 1}")
    axes[1, 0].plot(cutoff_ratios, sos[sec, 3, :], label=f"Section {sec + 1}")
    axes[1, 1].plot(cutoff_ratios, sos[sec, 4, :], label=f"Section {sec + 1}")
    axes[1, 2].plot(cutoff_ratios, sos[sec, 5, :], label=f"Section {sec + 1}")

axes[0, 0].set_ylabel("b0")
axes[0, 1].set_ylabel("b1")
axes[0, 2].set_ylabel("b2")
axes[1, 0].set_ylabel("a0")
axes[1, 1].set_ylabel("a1")
axes[1, 2].set_ylabel("a2")
for axes_row in axes:
    for ax in axes_row:
        ax.set_xscale("log")
        ax.legend()
        ax.grid(axis="both", which="major", lw=0.5, color="black")
        ax.grid(axis="x", which="minor", lw=0.2, color="black")
        ax.set_xlabel("Cutoff ratio $f_c/f_s$ [-]")

fig.suptitle(f"SOS coefs vs cutoff ratio\nfor order {order} Butterworth lowpass filter")
fig.tight_layout()
plt.show()
