import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint

K = 100000
m = 8
n = 5
ans = 0
for test in range(K):
  route = randint(0, 2, m + n)
  cnt = 0
  for i in range(m + n):
    if(route[i] == 0):
      cnt += 1
  if(cnt == m):
    ans += 1

print(ans / K)
