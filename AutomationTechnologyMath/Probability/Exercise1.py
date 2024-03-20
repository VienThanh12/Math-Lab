from re import T
N = 10000
case1 = 0
case2 = 0
for n in range(N): 
  result = randint(1, 7, 4)
  result2 = randint(1, 7, (2, 24))
  ok = False 
  for i in range(4):
    if(result[i] == 6):
      ok = True

  if(ok == True):
    case1 = case1 + 1
  ok = False
  for i in range(24):
    if result2[0][i] == 6 and result2[1][i] == 6:
      ok = True
  if(ok == True):
    case2 = case2 + 1
  
print("Simulation of Exercise 1")
print(case1 / N)
print(case2 / N)
