import numpy as np
import matplotlib.pyplot as plt
from numpy.random import permutation

K = 10000
ans = 0
ans_2 = 0
for test in range(K):
  s = permutation(52)
  np.random.shuffle(s)
  numbers=np.arange(1,14)
  numbers=np.hstack((numbers,numbers,numbers,numbers))
  suit=np.ones(13)
  suit=np.hstack((suit,2*suit,3*suit,4*suit))
  deck=np.zeros((52,2))
  deck[:,0]=numbers
  deck[:,1]=suit

  sdeck = deck[s,:]
  person1 = sdeck[:13]
  person2 = sdeck[13:26]
  person3 = sdeck[26:39]
  person4 = sdeck[39:52]
  ok1 = False
  ok2 = False
  ok3 = False
  ok4 = False
  for i in range(13):
    if(person1[i][0] == 1):
      ok1 = True
  for i in range(13):
    if(person2[i][0] == 1):
      ok2 = True
  for i in range(13):
    if(person3[i][0] == 1):
      ok3 = True
  for i in range(13):
    if(person4[i][0] == 1):
      ok4 = True 
  if ok1 == True and ok2 == True and ok3 == True and ok4 == True:
    ans += 1
  
  if(len(np.unique(person1[:, 1])) == 4 and len(np.unique(person2[:, 1])) == 4 and len(np.unique(person3[:, 1])) == 4 and len(np.unique(person4[:, 1])) == 4):
    ans_2 += 1
print("P(each person gets one ace (= number 1))")
print(ans / K)
print("P(each person gets cards from each 4 suites)")
print(ans_2 / K)




