import numpy as np
import matplotlib.pyplot as plt

K = 1

chosen_nums = 7
nums = 39

numbers = np.arange(1, nums, 1)

ans = 0

for test in range(K):
    chosen_numbers = permutation(numbers)[:chosen_nums]
    ok = True
    np.sort(chosen_numbers)

print(chosen_numbers)



