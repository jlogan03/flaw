# Calculate the expected gain an phase shift for testing SOS filters
import numpy as np
from scipy.signal import butter, freqz_sos

fc = 0.05  # cutoff frequency as a fraction of the sampling rate
sos = butter(4, fc, fs=1.0, btype="low", output="sos")
print(sos)
freqs, h = freqz_sos(sos, worN=4096, fs=1.0)
for fi in [0.01, 0.05, 0.1]:
    idx = np.argmin(np.abs(freqs - fi))
    gain = np.abs(h[idx])
    phase = np.angle(h[idx])
    print(f"Freq: {freqs[idx]:.4f}, Gain: {gain:.3e}, Phase: {phase:.3f} rad")
