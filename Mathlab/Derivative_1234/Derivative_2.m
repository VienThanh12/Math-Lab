clear
close all
L = 5
H = 3
B1 = -pi/9 % -20 degree
k1 =tan(B1)
A = [ L^2 , L;
	2 * L, 1]
B = [H;
k1]
sol = A\ B
a = sol(1)
b = sol(2)
delta = atan(b) * 180 / pi
g = 9.81
v0 = sqrt (-g / (2 * a)) / cos (delta * pi / 180)

f  = @(x) a .* x .^2 + b .*x

xmin = 0
xmax = L
x=xmin:(xmax-xmin)/100:xmax;
plot(x,f(x),'linewidth',2)

hold
grid
plot([L, L + 1], [H, H + k1],  'linewidth',2)
plot([5, 5], [0, 3] ,'linewidth',2)
plot([0, 5], [0, 0] ,'linewidth',2)
plot(5, 3)
title(['L = ', num2str(L), ', B =', num2str(-20), ': v0 = ', num2str(v0), ', a = ', num2str(delta)])
