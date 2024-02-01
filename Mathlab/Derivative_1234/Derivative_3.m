clear
close all
x1 = 1
y1 = 2
k1 = -0.25
m1 = 1
x2 = 5
y2 = 4
k2 = 0.5
m2 = -2


%coefficent matrix

A = [x1 ^ 5, x1 ^ 4, x1 ^ 3, x1 ^ 2, x1, 1;
       5 * (x1 ^ 4), 4 * (x1 ^ 3), 3 * (x1 ^ 2), 2 * x1, 1, 0;
       20 * (x1 ^ 3), 12 * (x1 ^ 2), 6 * x1, 2, 0, 0
       x2 ^ 5, x2 ^ 4, x2 ^ 3, x2 ^ 2, x2, 1;
       5 * (x2 ^ 4), 4 * (x2 ^ 3), 3 * (x2 ^ 2), 2 * x2, 1, 0;
       20 * (x2 ^ 3), 12 * (x2 ^ 2), 6 * x2, 2, 0, 0]
B = [y1;
       k1;
       m1;
       y2;
       k2;
       m2]
sol = A\ B
a = sol(1)
b = sol(2)
c = sol(3)
d = sol(4)
e = sol(5)
f = sol(6)

f=@(x) a * (x .^ 5) + b * (x .^ 4) + c  * (x .^ 3) + d * (x .^ 2) + e * x + f % f(x)
df=@(x) 5 * a * (x .^ 4) + 4 * b * (x .^ 3) + c * 3 * (x .^ 2) + 2 * d * x + e; %f'(x)
d2f=@(x) 20 * a * (x .^ 3) + 12 * b * (x .^ 2) + 6 * c * x + 2 * d; %f''(x)

k=@(x) d2f(x)./(sqrt(1+df(x).^2)).^3

x0=1
y0=f(x0)
df0=df(x0)
d2f0=d2f(x0)
k0=k(x0)

R0=abs(1/k0)

th0=atan2d(df0,1)

if k0>=0
    kx0=x0+R0*cosd(th0+90)
    ky0=y0+R0*sind(th0+90)
else
    kx0=x0+R0*cosd(th0-90)
    ky0=y0+R0*sind(th0-90)
end

t=0:1:360;
kx=kx0+R0*cosd(t);
ky=ky0+R0*sind(t);

% circle 2
x01=5
y01=4
df01=df(x01)
d2f01=d2f(x01)
k01=k(x01)

R01=abs(1/k01)

th01=atan2d(df01,1)

if k01>=0
    kx01=x01+R01*cosd(th01+90)
    ky01=y01+R01*sind(th01+90)
else
    kx01=x01+R01*cosd(th01-90)
    ky01=y01+R01*sind(th01-90)
end

t1=0:1:360;
kx1=kx01+R01*cosd(t1);
ky1=ky01+R01*sind(t1);

axis equal
% graph circle 1
figure(1)
hold on
xmin=1
xmax=5
x=xmin:(xmax-xmin)/100:xmax;
plot(x, f(x), 'linewidth', 1.5)

plot(kx,ky,'r','linewidth',1.5)
L=1 %horizontal width =2L

plot([x0,x0+L],[y0,y0+df0*L],'r','linewidth', 1.5)
plot([x0,kx0],[y0,ky0],'r','linewidth',1.5)
plot(x0,y0,'r.','markersize',20)
plot(kx0,ky0,'r.','markersize',10)

plot(kx1,ky1,'g','linewidth',1.5)

L=1 %horizontal width =2L

plot([x01,x01+L],[y01,y01+df01*L],'g','linewidth', 1.5)
plot([x01,kx01],[y01,ky01],'g','linewidth',1.5)
plot(x01,y01,'g.','markersize',20)
plot(kx01,ky01,'g.','markersize',10)
title(['k1', num2str(k1), ' , m1 = ', num2str(m1), ', k2 = ', num2str(k2), ', m2 = ', num2str(m2)])
hold off
grid


%graph 2
figure(2)
hold on
dx = (x2 - x1) / 100
x = x1:dx:x2;
plot(x1, k1, 'r.','markersize', 15)
plot(x2, k2, 'g.','markersize', 15)
dy = 5 * a * (x .^ 4) + 4 * b * (x .^ 3) + c * 3 * (x .^ 2) + 2 * d * x + e;
plot(x,dy,'linewidth', 1.5)
title(['dy(x) = 5 * a * (x .^ 4) + 4 * b * (x .^ 3) + c * 3 * (x .^ 2) + 2 * d * x + e'])
hold off
grid
%derivative
figure(3)
hold on
x = x1:dx:x2;
plot(x1, m1, 'r.','markersize', 15)
plot(x2, m2, 'g.','markersize', 15)
dyy = 20 * a * (x .^ 3) + 12 * b * (x .^ 2) + 6 * c * x + 2 * d;
plot(x, dyy, 'linewidth', 1.5)
title(['d2y(x) = 20 * a * (x .^ 3) + 12 * b * (x .^ 2) + 6 * c * x + 2 * d'])

hold off
grid


% kaarevuus
figure(4)
hold on
plot(x,k(x),'b','linewidth',1.5)
title('curvature \kappa(x)')
plot(x1, k(x1), 'r.','markersize', 15)
plot(x2, k(x2), 'g.','markersize', 15)
xlabel('x')
hold off
grid

