n=100000
X=np.random.normal(x1, sigma1,n) #n  N(mu,sigma^2)-distributed numbers

plt.plot(u, f, 'r', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x1,0,'ro')
plt.hist(X,bins=100,density=True,edgecolor='r',zorder=2)

plt.grid()
plt.show()
