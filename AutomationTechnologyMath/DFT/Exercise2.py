import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

U = 5
T = 3
N = 20

fs = N / T

f0 = fs / N
T_1 = 1 / f0
f = np.array([0, 1, 2, 3, 4, 5, 6])
time = np.linspace(0, T_1 / 2, 20000)
line = time * (U) / (T_1 / 2)

time_1 = np.linspace(T_1 / 2, T_1, 20000)
line_1 = U * 2 - time_1 * (U) / (T_1 / 2)

#########

n = np.arange(0, N / 2 + 1)
ts = n / fs
signal = ts * (U) / (T / 2)


n_2 = np.arange(N / 2 + 1, N)
ts_2 = n_2 / fs
signal_2 = U * 2 - ts_2 * (U) / (T / 2)

########
 
N_1 = N * 3 / 4

M_1 = int(N_1 / 2)

a_1 = np.zeros(M_1)
b_1 = np.zeros(M_1)

a0 = 1 / N_1 * np.sum(signal)

from numpy import cos, sin, pi
for k in range(1,M_1):
    a_1[k-1]=2/N_1*np.sum(signal*cos(2*pi*k*n/N_1))
    b_1[k-1]=2/N_1*np.sum(signal*sin(2*pi*k*n/N_1))

a_1[M_1-1]=1/N_1*np.sum(signal*cos(2*pi*M_1*n/N_1))
b_1[M_1-1]=1/N_1*np.sum(signal*sin(2*pi*M_1*n/N_1))

A_dft_1 = np.sqrt(a_1**2+b_1**2) #amplitudes

yt_1=0

for k in range(M_1):
    yt_1 = yt_1 + time * A_dft_1[k] / (T / 2)


####
 
N_2 = N * 3 / 4

M_2 = int(N_2 / 2)

a_2 = np.zeros(M_2)
b_2 = np.zeros(M_2)

from numpy import cos, sin, pi
for k in range(1,M_2):
    a_2[k-1]=2/N_2*np.sum(signal*cos(2*pi*k*n/N_2))
    b_2[k-1]=2/N_2*np.sum(signal*sin(2*pi*k*n/N_2))

a_2[M_2-1]=1/N_2*np.sum(signal*cos(2*pi*M_2*n/N_2))
b_2[M_2-1]=1/N_2*np.sum(signal*sin(2*pi*M_2*n/N_2))

A_dft_2 = np.sqrt(a_2**2+b_2**2) #amplitudes

yt_2=0

for k in range(M_1):
    yt_2 = yt_2 + time_1 * A_dft_2[k] / (T / 2)

yt_2 = 10 - yt_2


plt.subplot(311)

plt.plot(time, line, 'blue')
plt.plot(time_1, line_1, 'blue')
plt.grid()


plt.subplot(312)
plt.plot(ts, signal, 'ro')
plt.plot(ts_2, signal_2, 'ro')
plt.plot(time, yt_1,'blue', label = 'y(t)',lw=2)
plt.plot(time_1, yt_2,'blue',lw=2)
plt.legend()
# Plot the samples

plt.grid()

plt.subplot(313)
plt.bar(np.hstack((0,f)),np.hstack((a0,A_dft_1)),edgecolor='k',zorder=2)
plt.grid()
plt.xticks(np.hstack((0,f)))
plt.xlabel("frequency")

plt.tight_layout()

plt.ylabel("amplitudes")
plt.show()

