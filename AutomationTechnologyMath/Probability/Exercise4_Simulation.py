import numpy as np
import matplotlib.pyplot as plt
from numpy.random import permutation

K = 100000
n = 20


res = np.zeros(n + 1)

original_order = np.arange(1, n + 1)

for test in range(K):
  sp = np.arange(1, n + 1)
  ans = 0
  new_order = permutation(original_order)
  for i in range(n):
    if(original_order[i] == new_order[i]):
        ans += 1
  res[ans] += 1

for i in range(n):
  res[i] = res[i] / K

x=np.arange(0,n+1)
plt.figure(figsize=(8,4))
plt.bar(x,res,zorder=2)
plt.grid()
plt.xticks(x)
plt.xlabel('k')
plt.ylabel('P(k nroa pysyy paikoillaan)')
plt.title('n = ' + str(n))
plt.show()

