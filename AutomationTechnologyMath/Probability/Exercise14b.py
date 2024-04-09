n=100000
X=np.random.normal(x2, sigma2,n) #n  N(mu,sigma^2)-distributed numbers

plt.plot(u1, f1, 'g', label ='N($x_1, \sigma_1^2$)', linewidth = 2)
plt.plot(x2,0,'go')
plt.hist(X,bins=100,density=True,edgecolor='g',zorder=2)

plt.grid()
plt.show()
