n=100000
plt.subplot(311)
X=np.random.normal(x1, sigma1,n) #n  N(mu,sigma^2)-distributed numbers

plt.plot(u, f, 'r', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x1,0,'ro')
plt.hist(X,bins=100,density=True,edgecolor='k',zorder=2, facecolor = 'r')

plt.grid()

plt.subplot(312)
n=100000
X=np.random.normal(x2, sigma2,n) #n  N(mu,sigma^2)-distributed numbers

plt.plot(u1, f1, 'g', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x2,0,'go')
plt.hist(X,bins=100,density=True,edgecolor='k',zorder=2, facecolor = 'g')

plt.grid()

plt.subplot(313)
n=100000
X=np.random.normal(x, sigma,n) #n  N(mu,sigma^2)-distributed numbers

plt.plot(u2, f2, 'b', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x,0,'ro')
plt.hist(X,bins=100,density=True,edgecolor='k',zorder=2, facecolor = 'b')

plt.grid()
plt.show()
