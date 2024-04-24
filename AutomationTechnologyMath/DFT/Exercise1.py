import numpy as np
import matplotlib.pyplot as plt

# Define parameters
A = 3  # Amplitude
f = 22  # Frequency in Hz
phi = 0 # Phase angle in radians
N = 20  # Number of samples
fs = 100  # Sampling frequency in Hz

n = np.arange(0, N)
ts = n / fs

f0 = fs / N
T = 1 / f0
t = np.linspace(0, T, 200000)

signal = A * np.cos(2 * pi * f * ts + phi)

t_2 = np.linspace(0, T, 20000)
signal_2 = A * np.cos(2 * np.pi * f * t_2 + phi)

print(signal)

######### 

M = int(N / 2)

a0 = 1 / N * np.sum(signal)
a = np.zeros(M)
b = np.zeros(M)

from numpy import cos, sin, pi
for k in range(1,M):
    a[k-1]=2/N*np.sum(signal*cos(2*pi*k*n/N))
    b[k-1]=2/N*np.sum(signal*sin(2*pi*k*n/N))

a[M-1]=1/N*np.sum(signal*cos(2*pi*M*n/N))
b[M-1]=1/N*np.sum(signal*sin(2*pi*M*n/N))

A_dft = np.sqrt(a**2+b**2) #amplitudes
phi_dft = np.arctan2(-b,a) #phase angles

f_dft=np.arange(f0,(M+1)*f0,f0) #frequencies

plt.figure(figsize=(7,4))

yt=a0
for k in range(M):
    yt=yt+A_dft[k]*cos(2*pi*f_dft[k]*t+phi_dft[k])

# Plot the signal f(t)
plt.subplot(211)
plt.plot(t, yt,'orange', label = 'y(t)',lw=2)

plt.plot(t_2, signal_2, label='f(t)')

plt.legend(fontsize = 14)
plt.xlabel('time',fontsize=14)
plt.title('A = ' + str(A) + ', f = '+ str(f) + ', $\phi$ = '+str(phi) + ', N = ' + str(N)+ ', $f_s$ = '+ str(fs))

plt.plot(ts, signal, 'ro')

# Plot the samples

plt.grid(True)

plt.tight_layout()

plt.subplot(212)
plt.bar(np.hstack((0,f_dft)),np.hstack((a0,A_dft)),edgecolor='k',zorder=2)
plt.grid()
plt.xticks(np.hstack((0,f_dft)))
plt.xlabel("frequency")
plt.ylabel("amplitudes")
plt.show()

plt.show()
