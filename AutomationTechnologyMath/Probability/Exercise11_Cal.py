import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint
from scipy.special import comb

n = 5

res = np.zeros(10)
ans = 0
p = 0.4

for m in range(1, 11):
  ans = 0
  for r in range(m, m + n, 1):
    ans = ans + comb(m + n - 1, r) * (p ** r) * ((1 - p) ** (m + n - 1 - r))
  res[m - 1] = ans

res_1 = np.zeros(10)
ans = 0
p = 0.5

for m in range(1, 11):
  ans = 0
  for r in range(m, m + n, 1):
    ans = ans + comb(m + n - 1, r) * (p ** r) * ((1 - p) ** (m + n - 1 - r))
  res_1[m - 1] = ans

res_2 = np.zeros(10)
ans = 0
p = 0.6

for m in range(1, 11):
  ans = 0
  for r in range(m, m + n, 1):
    ans = ans + comb(m + n - 1, r) * (p ** r) * ((1 - p) ** (m + n - 1 - r))
  res_2[m - 1] = ans

x = np.arange(1, 11, 1)
plt.figure(figsize = (8,5))
plt.plot(x, res_1, 'g.-', label = 'p = 0.5' ,markersize = 12)
plt.plot(x, res, 'r.-', label = 'p = 0.4',markersize = 12)
plt.plot(x, res_2, 'b.-', label = 'p = 0.6' ,markersize = 12)
plt.grid()
plt.xlabel('m', fontsize = 15)
plt.ylabel('P_mn', fontsize = 15)
plt.legend(fontsize = 15)
plt.title('n = '+str(n), fontsize = 15)
plt.show()
