import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint

K = 100000
n = 20

graph = np.zeros(n)
for test in range(K):
  ans = 0
  res = 0
  temp = 0
  for i in range(20):
    coin = randint(1, 3)
    if(temp == coin):
      res = res + 1
      ans = max(res, ans)
    else:
      res = 0
    temp = coin
  graph[ans] += 1


x = np.arange(0, n)
plt.figure(figsize = (6, 4))
plt.bar(x, graph / K, zorder = 1)
plt.grid()
plt.xlabel('k')
plt.ylabel('P(k consecutive heads)')
plt.title('n = ' + str(n))
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.show()
print(p)
