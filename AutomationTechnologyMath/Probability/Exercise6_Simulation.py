# Simulation
import numpy as np
import matplotlib.pyplot as plt

K = 100000

choose_num = 10
win_num = 20
non_num = 50

ans = np.zeros(choose_num)
numbers=np.arange(1, win_num + non_num, 1)

winning_numbers = permutation(numbers)[:win_num]

for test in range(K):
  choose_numbers = permutation(numbers)[:choose_num]
  winning = np.intersect1d(choose_numbers,winning_numbers)
  k = winning.size
  ans[k] += 1

plt.figure(figsize = (8, 5))
x = np.arange(0, 10, 1)
plt.bar(x, ans / K, label = 'simulointi', zorder = 2)
plt.plot(x, p, 'r.',label = 'tarkka', markersize = 15)
plt.legend(fontsize=14)
plt.xlabel('k')
plt.ylabel('k oikein kierrosten osuus')
plt.title('kierroksia K = '+str(K))
plt.grid()
plt.show()
