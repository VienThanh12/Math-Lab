import matplotlib.pyplot as plt
from sklearn import datasets
import numpy as np
#%%
#X,y=datasets.load_wine(return_X_y=True)

#X = 178x13-matrix, datapoints as rows
#y=0,1,2=viini Barolo,Grignolino,Barbera
#%%
X,y=datasets.fetch_olivetti_faces(return_X_y=True)

#X=400x4096-matrix, datapoints as rows
#y=0,1,2,...,39=person number

plt.figure(figsize=(8,8))

for i in range(40):
  plt.subplot(5, 8, i + 1)
  plt.imshow(X[i * 10, :].reshape(64,64), cmap='gray')
  plt.title('y = {}'.format(i))
  plt.axis('off')

#scaling
mu = np.mean(X, axis = 0)
Xhat = X - mu

U,s,V = np.linalg.svd(Xhat,full_matrices=False)
Z = U * s

plt.figure(figsize=(8,4))
plt.plot(s)
plt.grid()
plt.xlabel('k',fontsize=14)
plt.ylabel('singulaariarvo $\sigma_k$',fontsize=14)
plt.show()
u = np.cumsum(s**2)
v = np.sum(s**2)
plt.figure(figsize=(8,4))
plt.plot(u / v)
plt.grid()
plt.xlabel('k', fontsize=14)
plt.title('varience of the principal component 1...k', fontsize=14)
plt.show()

plt.subplot(5, 8, 1)
plt.imshow(mu.reshape(64,64), cmap='gray')
plt.axis('off')
plt.title('u')

for i in range(1, 32):
  plt.subplot(5, 8, i + 1)
  plt.imshow(V[i, :].reshape(64,64), cmap='gray')
  plt.title('v = {}'.format(i))
  plt.axis('off')

Z = U * s
plt.figure(figsize=(8, 8))

plt.subplot(3, 7, 1)
la = mu + Z[65, :20].reshape(1, 20)@V[:20,:]
plt.imshow(la.reshape(64, 64), cmap='gray')
plt.title('k = {}'.format(6))
plt.axis('off')
i = 2
for k in range(20, 420, 20):
  plt.subplot(3, 7, i)
  la = mu + Z[65,:k].reshape(1,k)@V[:k,:]
  plt.imshow(la.reshape(64, 64), cmap='gray')
  plt.title('k = {}'.format(k))
  i += 1
  plt.axis('off')
