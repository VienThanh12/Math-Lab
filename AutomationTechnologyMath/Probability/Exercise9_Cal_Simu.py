import numpy as np
import matplotlib.pyplot as plt
from numpy.random import permutation
from scipy.special import comb
from numpy.random import randint

K = 10000
questions = 20

# ((20) C k) * (3 ^ (20 - k))
# Ptotal = 4 ^ 20


P_total = 4 ** (questions)


p = np.zeros(questions) #p[k] = P(k winning numbers), k=0,1,2,...,7
#from scipy.special import comb
for k in range(questions):
  p[k] = comb(questions, k) * (3 ** (questions - k)) / P_total

res = np.zeros(questions)

for test in range(K):
  ans = 0
  for i in range(questions):
    x = randint(1, 5)
    y = randint(1, 5)
    if(x == y):
      ans += 1
  res[ans] += 1

x = np.arange(0, questions)
plt.figure(figsize = (6, 4))
plt.bar(x, p, label = 'tarkka', zorder = 1)
plt.plot(x, res / K, 'r.',label = 'simu', markersize = 15)
plt.grid()
plt.xlabel('k')
plt.ylabel('P(k oikein)')
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)
plt.legend(fontsize = 14)
plt.show()
print(p)
