import numpy as np
import matplotlib.pyplot as plt

K = 100000

n = 10
t = 50

ans = 0
for k in range(K):
  check = np.zeros(n + 1)
  coupons = randint(1, n + 1, t)
  for i in range(t):
    check[coupons[i]] += 1
  ok = True
  for i in range(1, n + 1):
    if(check[i] == 0):
      ok = False
  if(ok == True):
    ans += 1
  
print(ans / K)
