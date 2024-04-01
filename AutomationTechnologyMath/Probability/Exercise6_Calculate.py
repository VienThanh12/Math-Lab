from scipy.special import comb
import numpy as np
import matplotlib.pyplot as plt

choose_num = 10
win_num = 20
non_num = 50

# player chooses 10 num from 1 -> 70 (we have 50 non win and 20 win)
# So the total of posibilities will be win_num + non_num C choosenum

P_total = comb(win_num + non_num, choose_num)

p=np.zeros(choose_num) #p[k] = P(k winning numbers), k=0,1,2,...,7
#from scipy.special import comb
for k in range(choose_num):
  p[k] = comb(win_num, k) * comb(non_num, choose_num - k) / P_total

x = np.arange(0,10)
plt.figure(figsize = (6, 4))
plt.bar(x, p, label = 'simulointi', zorder = 2)
plt.grid()
plt.xlabel('k')
plt.ylabel('P(k winning numbers)')
plt.plot(x, p,'r.', label = 'tarkka',markersize = 15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
plt.show()
print(p)

