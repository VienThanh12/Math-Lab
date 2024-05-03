import matplotlib.pyplot as plt
from sklearn import datasets
import numpy as np
#%%
X,y = datasets.load_wine(return_X_y = True)

#X = 178x13-matrix, datapoints as rows
#y=0,1,2=viini Barolo,Grignolino,Barbera
#%%
#X,y=datasets.fetch_olivetti_faces(return_X_y=True)

#X=400x4096-matrix, datapoints as rows
#y=0,1,2,...,39=person number

#scaling
mu = np.mean(X, axis = 0)
stds = np.std(X, axis = 0)

Xhat = (X - mu) / stds

U,s,V = np.linalg.svd(Xhat, full_matrices = False)
Z = U * s

'''
for k in range(4):
    plt.subplot(1,4,k+1)
    plt.hist(Xhat[:,k],zorder=2,edgecolor='k')
    plt.grid()
    plt.title('x'+str(k+1),fontsize=14)
'''

plt.figure(figsize = (8,4))
plt.plot(s)
for i in range(13):
  plt.plot(i, s[i], 'co')
plt.grid()
plt.xlabel('k',fontsize=14)
plt.ylabel('singulaariarvo $\sigma_k$',fontsize=14)
plt.show()

u = np.cumsum(s ** 2)
v = np.sum(s ** 2)
plt.figure(figsize = (8,4))

z = u / v
plt.plot(u / v)
plt.grid()
for i in range(13):
  plt.plot(i, z[i], 'co')
plt.xlabel('k', fontsize = 14)
plt.title('varience of the principal component 1...k', fontsize=14)
plt.show()

plt.plot(Z[y==0,0],Z[y==0,1],'r.',label='y=0')
plt.plot(Z[y==1,0],Z[y==1,1],'g.',label='y=1')
plt.plot(Z[y==2,0],Z[y==2,1],'b.',label='y=2')
plt.grid()
plt.axis('equal')
plt.xlabel('z1')
plt.ylabel('z2',rotation=0)
plt.legend()
plt.title('picture using principal components z1 and z2')
plt.show()
