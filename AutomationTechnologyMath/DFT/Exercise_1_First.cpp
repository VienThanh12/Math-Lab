import numpy as np
import matplotlib.pyplot as plt

# Define parameters
A = 3  # Amplitude
f = 20  # Frequency in Hz
phi = 0 # Phase angle in radians
N = 20  # Number of samples
fs = 100  # Sampling frequency in Hz

n = np.arange(0, N)
ts = n / fs

f0 = fs / N
T = 1 / f0
t = np.linspace(0, T, 20)

signal = A * np.cos(2 * pi * f * ts + phi)

t_2 = np.linspace(0, T, 20000)
signal_2 = A * np.cos(2 * np.pi * f * t_2 + phi)

print(signal)

# Plot the signal f(t)
plt.subplot(2, 1, 1)
plt.plot(t_2, signal_2)
plt.title('Original Signal')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.plot(ts, signal, 'ro')

# Plot the samples

plt.grid(True)

plt.tight_layout()
plt.show()
