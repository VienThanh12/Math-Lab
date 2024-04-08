import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint
from scipy.special import comb

#GIVEN
x1 = 1
sigma1 = 1.5
x2 = 5
sigma2 = 2.5


# Solve

t = sigma1 ** 2 + sigma2 ** 2
a1 = sigma2 * sigma2 / t
a2 = sigma1 * sigma1 / t
x = a1 * x1 + a2 * x2

sigma = np.sqrt(sigma1 ** 2 * sigma2 ** 2 / t)

u=np.linspace(x1-4*sigma1,x1+4*sigma1,100)
f=1/(np.sqrt(2*np.pi)*sigma1)*np.exp(-1/2*((u-x1)/sigma1)**2)


u1=np.linspace(x2-4*sigma2,x2+4*sigma2,100)
f1=1/(np.sqrt(2*np.pi)*sigma2)*np.exp(-1/2*((u1-x2)/sigma2)**2)

u2=np.linspace(x-4*sigma,x+4*sigma,100)
f2=1/(np.sqrt(2*np.pi)*sigma)*np.exp(-1/2*((u2-x)/sigma)**2)


plt.figure(figsize=(8,5))
plt.plot(u, f, 'r', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x1,0,'ro')

plt.plot(u1, f1, 'g', label ='N($x_2, \sigma_2^2$)', linewidth = 2)
plt.plot(x2,0,'go')

plt.plot(u2, f2, 'b', label ='N($x, \sigma^2$)', linewidth = 2)
plt.plot(x,0,'bo')


plt.grid()
plt.legend(fontsize = 14)
#plt.title('x_1 = '+str(x1)+', $\sigma_1$ = '+str(sigma1), 'x_2 = '+str(x2)+', $\sigma_2$ = '+str(sigma2),
#          fontsize=14)
plt.legend(fontsize=14)
plt.show()
