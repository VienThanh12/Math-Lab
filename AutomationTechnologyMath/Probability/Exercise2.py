#%% simulation
K = 100000#simulation rounds
n = 100 #number of persons

ans = 0

for k in range(K):
  sp = np.random.randint(1, 366, n)
  check = np.zeros(367)
  for i in range(n):
    check[sp[i]] += 1
  res = 0
  for i in range(n):
    res = max(res, check[sp[i]])
  if(res >= 3):
    ans += 1

print(ans / K)
