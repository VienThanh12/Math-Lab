import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint
from scipy.special import comb

K = 100000
ans = 0
lose = 0
for test in range(K):
  sum = 0
  x = randint(1, 7, 1)
  y = randint(1, 7, 1)
  sum = x + y
  if(sum == 7 or sum == 11):
    ans += 1
  elif(sum == 2 or sum == 3 or sum == 12):
    lose += 1
  else: 
    res = 0 
    while(res != sum and res != 7):
      x = randint(1, 7, 1)
      y = randint(1, 7, 1)
      res = x + y
    if(res == sum):
      ans += 1
print(ans / K)
