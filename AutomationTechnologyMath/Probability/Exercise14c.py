n=100000
X=np.random.normal(x, sigma,n) #n  N(mu,sigma^2)-distributed numbers

plt.plot(u2, f2, 'b', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x,0,'ro')
plt.hist(X,bins=100,density=True,edgecolor='b',zorder=2)

plt.grid()
plt.show()
