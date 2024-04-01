from scipy.special import comb

m = 5
n = 7
# Choose m from m + n so we will have m+nCm
# Each of the turn, the traveller moves up or right. So it will be 2 ^ (m + n) case
P_total = 2 ** (m + n)

print(comb(m + n, n) / (P_total))
