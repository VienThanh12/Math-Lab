import numpy as np
import matplotlib.pyplot as plt

# GIVEN

mu_x1 = 2
mu_y1 = 3
mu_1 = np.array([[mu_x1, mu_y1]]).T

mu_x2 = 1.5
mu_y2 = 2
mu_2 = np.array([[mu_x2, mu_y2]]).T

sigma_x1 = 1
sigma_y1 = 1.5

p_1 = 0.3
sigma_xy1 = p_1 * sigma_x1 * sigma_y1

sigma_x2 = 1.5
sigma_y2 = 0.8

p_2 = 0.4

sigma_xy2 = p_2 * sigma_x2 * sigma_y2

# Solve

def ellipse2D(mu, S, p):
    vals, vect = np.linalg.eig(S)
    la1 = vals[0]
    la2 = vals[1]
    u = vect[:,0]
    v = vect[:,1]
    a = np.sqrt(p*la1)
    b = np.sqrt(p*la2)
    t = np.linspace(0,2*np.pi,100)
    co = np.cos(t)
    si = np.sin(t)
    x = mu[0,0] + a * co * u[0] + b * si * v[0]
    y = mu[1,0] + a * co * u[1] + b * si * v[1]
    return x,y

S_1 = np.array([[sigma_x1 ** 2, sigma_xy1],
            [sigma_xy1, sigma_y1 ** 2]])

S_2 = np.array([[sigma_x2 ** 2, sigma_xy2],
            [sigma_xy2, sigma_y2 ** 2]])

A_1 = S_2@np.linalg.inv(S_1 + S_2)
A_2 = S_1@np.linalg.inv(S_1 + S_2)
mu_3 = A_1@mu_1 + A_2@mu_2
S_3 = S_2@np.linalg.inv(S_1 + S_2)@S_1
sigma_x3 = np.sqrt(S_3[0, 0])
sigma_y3 = np.sqrt(S_3[1, 1])

p = 6.0
x1,y1 = ellipse2D(mu_1, S_1, p) 
x2,y2 = ellipse2D(mu_2, S_2, p) 
x3,y3 = ellipse2D(mu_3, S_3, p)


plt.plot(x1, y1, 'r', label = 'N(X1)')
plt.plot(mu_x1, mu_y1, 'ro', markersize = 10)
plt.plot(x2,y2,'g',label = 'N(X2)')
plt.plot(mu_x2, mu_y2, 'go', markersize = 10)
plt.plot(x3, y3, 'b', label = 'N(X3)')
plt.plot(mu_3[0], mu_3[1], 'bo', markersize = 10)
plt.grid()
plt.legend(fontsize = 15)

plt.title('$\sigma_x1$ = '+str(sigma_x1)+', $\sigma_y1$ = '+str(sigma_y1)+ r', $\rho1$ = '+str(p_1)+
          ', $\sigma_x2$ = '+str(sigma_x2)+', $\sigma_y2$ = '+str(sigma_y2)+ r', $\rho2$ = '+str(p_2),fontsize=14)
plt.show()


N = 10000
x_1m = np.random.normal(mu_x1, sigma_x1, N)
y_1m = np.random.normal(mu_y1, sigma_y1, N)

plt.plot(x_1m, y_1m, 'r.', markersize = 1)
plt.plot(x1, y1, 'r', label='N(X1)')
plt.plot(mu_x1, mu_y1, 'ro', markersize = 5, label='')
plt.grid()
plt.legend(fontsize = 15)
plt.show()

N = 10000
x_2m = np.random.normal(mu_x2, sigma_x2, N)
y_2m = np.random.normal(mu_y2, sigma_y2, N)

plt.plot(x_2m, y_2m, 'g.', markersize = 1)
plt.plot(x2, y2, 'g', label='N(X2)')
plt.plot(mu_x2, mu_y2, 'go', markersize = 5, label='')
plt.grid()
plt.legend(fontsize = 15)
plt.show()
     

N = 10000
x_m = np.random.normal(mu_3[0], sigma_x3, N)
y_m = np.random.normal(mu_3[1], sigma_y3, N)

plt.figure(figsize = (5,5))
plt.plot(x_m, y_m, 'b.', markersize = 1)
plt.plot(x3, y3, 'b', label='N(X3)')
plt.plot(mu_3[0],mu_3[1],'bo',markersize = 5,label='')
plt.grid()
plt.legend(fontsize = 15)
plt.show()
