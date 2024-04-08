import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint
from scipy.special import comb

K = 10000
ans = 0

for test in range(K):
  x = randint(1, 60, 1)
  y = randint(1, 60, 1)
  if(x > y - 10 and x < y + 10):
    ans += 1
print(ans / K)
