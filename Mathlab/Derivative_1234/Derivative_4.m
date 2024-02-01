clear
close all
A = 1
a = 1.5
u = 1.57

B = 1.5
b = 1
T = 12.57

x=@(t) A .* sin(a .* t + u)
y=@(t) B*sin(b .* t)

dx=@(t) A .* a .* cos(a * t + u)
dy=@(t) B .* b .* cos(b * t)

d2x= @(t) -A .* a .^ 2 * sin(a .* t + u)
d2y= @(t) -B .* b .^ 2 * sin(b .* t)
% value

t0 = 4.5

% Curvature

k = @(t) (dx(t) .* d2y(t) - d2x(t) .* dy(t)) ./ ((sqrt(dx(t) .^ 2 + dy(t) .^ 2)) .^ 3)

% value

x0 = x(t0)
y0 = y(t0)
dx0 = dx(t0)
dy0 = dy(t0)
d2x0 = d2x(t0)
d2y0 = d2y(t0)
k0 = k(t0)

% radius
R0 = 1 / abs(k0)
th0=atan2d(dy0, dx0)

%center of the circle of curvature
if k0>=0
    kx0=x0+R0*cosd(th0+90)
    ky0=y0+R0*sind(th0+90)
else
    kx0=x0+R0*cosd(th0-90)
    ky0=y0+R0*sind(th0-90)
end
%points on the circle of curvature
t=0:1:360;
kx=kx0+R0*cosd(t);
ky=ky0+R0*sind(t);



axis equal
t = 0:T/100:T;
figure(1)
hold on
plot(x(t),y(t),'b','linewidth',1.5)
plot(kx,ky,'g','linewidth',1.5)

L=0.5 %horizontal width =2L
plot([x0,kx0],[y0,ky0],'r','linewidth',1.5)
plot(x0,y0,'r.','markersize',10)
plot(kx0,ky0,'g.','markersize',10)

plot([x0,x0+dx0],[y0,y0+dy0],'r','linewidth',1.5)
plot(x0,y0,'r.','markersize',15)
title(['A = ', num2str(1), ', a = ', num2str(a), ', u = ', num2str(u), ', B = ', num2str(B), ', b = ', num2str(b), ', T = ', num2str(T)
 ,'t = ' , num2str(t0), 'P(t) = [', num2str(x0), ',', num2str(y0), ']','T(t) = [', num2str(dx0), ' , ', num2str(dy0), ']', 'k(t) = ', num2str(k0),' R(t) = ', num2str(R0)
])

grid
hold off
% length of the tangent vector |T(t)|
axis equal
TT = @(t) sqrt(dx(t) .^ 2 + dy(t) .^ 2)
figure(2)
grid
hold on
t = 0:T/100:T;
plot(t / pi, TT(t), 'linewidth', 1.5)
title('length of the tangent vector |T(t)|')
xlabel('t / π')
hold off
% curvature


hold on
figure(3)

grid
title('f(t)')
t = 0:T/100:T;
plot(t / pi, k(t),'b','linewidth',1.5)
grid
title('curvature \kappa(t)')
xlabel('t / π')





