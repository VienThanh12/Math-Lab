import numpy as np
import matplotlib.pyplot as plt
from math import factorial
from numpy.random import randint

result1 = randint(1,7,4)

P1_none = 5 ** 4 / 6 ** 4
P1 = 1 - P1_none

P2_none = (35) ** 24 / 36 ** 24
P2 = 1 - P2_none
print("Exercise 1: Calculate the probabilites")
print(P1)
print(P2)
