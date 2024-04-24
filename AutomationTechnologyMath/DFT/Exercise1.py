import numpy as np
import matplotlib.pyplot as plt

#Given

A_ft = 3
# First case
#f_ft = 22
# Second case
f_ft = 20
phi_ft = 0
N_ft = 20

fs= 100 #sampling frequency

N = 20
#Solve
t=np.linspace(0,T,1000)

#First case
# y = np.array([3, 0.5, -2.8, -1.5, 2.2, 2.4, -1.2, -2.9, 0.2, 3, 0.9, -2.7, -1.9, 1.9, 2.6, -0.9, -3, -0.2, 2.9, 1.3])
# Second case
y = np.array([3, 0.9, -2.4, -2.4, 0.9, 3, 0.9, -2.4, -2.4, 0.9, 3, 0.9, -2.4, -2.4, 0.9, 3, 0.9, -2.4, -2.4, 0.9])
M = int(N / 2)
n = np.arange(0, N)

a0 = 1 / N * np.sum(y)
a = np.zeros(M)
b = np.zeros(M)

from numpy import cos, sin, pi
for k in range(1,M):
    a[k-1]=2/N*np.sum(y*cos(2*pi*k*n/N))
    b[k-1]=2/N*np.sum(y*sin(2*pi*k*n/N))

a[M-1]=1/N*np.sum(y*cos(2*pi*M*n/N))
b[M-1]=1/N*np.sum(y*sin(2*pi*M*n/N))

A=np.sqrt(a**2+b**2) #amplitudes
phi=np.arctan2(-b,a) #phase angles

f0=fs/N #basic frequency
f=np.arange(f0,(M+1)*f0,f0) #frequencies

T=1/f0
ts=n/fs
print(ts)

t=np.linspace(0,T,1000)

plt.figure(figsize=(7,4))

yt=a0
for k in range(M):
    yt=yt+A[k]*cos(2*pi*f[k]*t+phi[k])

plt.subplot(211)
plt.plot(t, A_ft * cos(2 * pi * f_ft * t + phi_ft),lw=2,label='f(t)')
plt.plot(t, yt,'orange', label = 'y(t)',lw=2)
plt.plot(ts,y,'ro')
plt.grid()
plt.legend(fontsize = 14)
plt.xlabel('time',fontsize=14)
plt.title('A = ' + str(A_ft) + ', f = '+ str(f_ft) + ', $\phi$ = '+str(phi_ft) + ', N = ' + str(N_ft)+ ', $f_s$ = '+ str(fs))

plt.subplot(212)
plt.bar(np.hstack((0,f)),np.hstack((a0,A)),edgecolor='k',zorder=2)
plt.grid()
plt.xticks(np.hstack((0,f)))
plt.xlabel("frequency")
plt.ylabel("amplitudes")
plt.show()
