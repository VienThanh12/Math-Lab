% Ha Vien (k2202863)
% Vien.Ha@seamk.fi
% Differential And Integral calculus
% 7/1/2024

% LIST OF THE EXERCISES
% Derivate: Exercise 1, 2, 3, 4 (4/4)
% Extrenum: Exercise 1, 2, 3, 4, 5 (5/5)
% Position Velocity Acceleration: Exercise 1, 2, 3 (3/7)
% Acceleration Velocity Position: Exercise 1, 2, 3, 4 (4/4)
% Spline and Bézier curves: Exercise 1 (1/4)
% Integral: Exercise 1, 2, 3, 4, 5, 6, 7 (7/7)
% Taylor Polynomial: Exercise 1 (1/1)

% TOTAL EXERCISES I HAVE COMPLEDTED: 25 / 32

% DERIVATIVE
% Exercise 1

clear
close all

f =@(x) 0.3.*x.^3 - 0.5.*x.^2 - 3.*x + 2
df=@(x) 0.9.*x.^2-x-3

x0=1.6
f0=f(x0)
df0=df(x0)
dx=f(x0)/df(x0)
x1=x0-dx
x=0:5/100:5;
plot(x,f(x),'linewidth',1.5)
hold
%tangent
L=dx %width = 2*L
plot([x0-L,x0+L],[f(x0)-df(x0)*L,f(x0)+df(x0)*L],'k','linewidth',1.5)
plot(x0,f(x0),'r.','markersize',20)
plot(x1,0,'g.','markersize',20)
grid
title(['x0 = ', num2str(x0), ', x1 = ', num2str(x1)])
hold off

% Exercise 2

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

% Exercise 3

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

% Exercise 4

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

% EXTRENUM
% Exercise 1

clear
close all

E = 1
t = 15
k = 3
f = @(R) ((E .* R) / (R + k) - (E .* R) / (R + t));

% df = -(E*(t-k)*(R^2-k*t))/((R+k)^2*(R+t)^2)
% R = - sqrt(k * t) or sqrt(k * t)
% R > 0
% R = sqrt(k * t sub to the equation
R = sqrt(k *t) % R max
f(sqrt(k * t)) % f max

x=0:0.01:60;

u = @(x) ((E .* x) ./ (x + k) - ((E .* x) ./ (x + t)))

%[umax,indmax]=max(u(x)) %largest value of vector f is its element number indmin
%xmax=x(indmax)

hold on
plot(x, u(x), 'linewidth',1.5)
grid
plot(R,f(R),'r.','markersize',15);
title(['E = ', num2str(E), ', Rt = ', num2str(t), ', Rk = ', num2str(k), ' , max u(', num2str(R), ') = ', num2str(f(R))])
ylabel('jännite u')
xlabel('R')
hold off

% Exercise 2

clear
close all

v = 50
e = 100
r = 10000

c = @(q) v .* (q ./ 2) + e .* (r ./ q)

% dc = v/2-(e*r)/q^2
% q = sqrt(2)*sqrt(e)*sqrt(r/v) or -sqrt(2)*sqrt(e)*sqrt(r/v)
% q > 0

u = sqrt(2)*sqrt(e)*sqrt(r/v)
c(u)
hold on
q = 0:10:800

%[cmin,indmin] = min(c(q)) %min value of vector f is its element number indmin
%qmin = q(indmin)

plot(q, c(q), 'linewidth', 1.5)
plot(u,c(u),'r.','markersize',15);
xlabel('eräkoko q')
title(['r = ', num2str(r), ', kv = ', num2str(v), ', ke = ', num2str(e), ': min = ', num2str(c(u)), ', kun q = ', num2str(u)])
ylabel('vuotuiset kustammukset')
grid
hold off

% Exercise 3

clear
close all

% MA + MB = 2 * sqrt(h ^ 2 + (b/2) ^ 2)
% MC = sqrt(a ^ 2 - (b/2) ^ 2) - h
% MA + MB + MC = 2 * sqrt(h ^ 2 + (b/2) ^ 2) + sqrt(a ^ 2 - (b/2) ^ 2) - h


a = 5
b = 6

f = @(h) 2 .* sqrt(h .^ 2 + (b ./ 2) .^ 2) + sqrt(a .^ 2 - (b ./ 2) .^ 2) - h
% df = (2*h)/sqrt(h^2+b^2/4)-1
% df = 0 => h = sqrt(b ^ 2 / 12)
u = sqrt(b ^ 2 / 12)
f(u)

k = sqrt(a ^ 2 - (b / 2) ^ 2)

figure(1)

hold on
h = 0:0.01:k


%[fmin,indmin] = min(f(h)) %min value of vector f is its element number indmin
%hmin = h(indmin)

if u >= k
    u = k
end

plot(h, f(h), 'linewidth', 1.5)
plot(u, f(u),'r.','markersize',15)
xlabel('h')
ylabel('s')
title(['a = ', num2str(a), ', b = ', num2str(b), ', h0 = ', num2str(u)])
grid
hold off

figure(2)
hold on

plot([-b / 2, b / 2], [0, 0], 'b','linewidth',1.5)

if( u != k)
    plot([-b / 2, 0], [0, sqrt(a ^ 2 - (b / 2) ^ 2)], 'b','linewidth',1.5)
    plot([b / 2, 0], [0, sqrt(a ^ 2 - (b / 2) ^ 2)], 'b','linewidth',1.5)
    plot(0, sqrt(a ^ 2 - (b / 2) ^ 2), 'b.' ,'markersize', 15)
end

plot(0, u, 'r.' ,'markersize', 17)

plot([0, 0], [u, sqrt(a ^ 2 - (b / 2) ^ 2)], 'r','linewidth',1.7)
plot([0, -b/2], [u, 0], 'r','linewidth',1.7)
plot([0, b/2], [u, 0], 'r','linewidth',1.7)
plot(-b/2, 0, 'b.' ,'markersize', 15)
plot(b/2, 0, 'b.' ,'markersize', 15)

grid
hold off

% Exercise 4

clear
close all

% PQ = sqrt(h ^ 2 + x ^ 2)
% QB = L - x

% t(x)  = sqrt(h ^ 2 + x ^ 2) + (L - x) / v

L = 10
h = 4
v = 1.05

t = @(x) sqrt(h .^ 2 + x .^ 2) + (L - x) ./ v
% t' = x/sqrt(x^2+h^2)-1/v
% t' = 0
% x = sqrt(h ^ 2 / (v ^ 2 - 1))
w = sqrt(h ^ 2 / (v ^ 2 - 1))
t(w)

if (w > L)
  w = L
end

figure(1)
hold on

x = 0:0.01:L

%[tmin,indmin] = min(t(x)) %min value of vector f is its element number indmin
%xmin = x(indmin)

plot(x, t(x), 'linewidth', 1.5)
plot(w, t(w),'r.','markersize',15)
ylabel('aika t reitillä PQB')
xlabel('AQ')
title(['L = ', num2str(L), ', h = ', num2str(h), ', v = ', num2str(v), ', min t(', num2str(w), ') = ', num2str(t(w))])
grid
hold off

figure(2)
hold on
plot(0, h, 'b.' ,'markersize', 15)
plot(w, 0, 'b.' ,'markersize', 15)

plot([0, w], [h, 0], 'b','linewidth',1.7)
plot([L, w], [0, 0], 'b','linewidth',1.7)
plot([0, w], [0, 0], 'black','linewidth',1)
title(['nopein retti, AQ = ', num2str(w)])
ylim([-1,5])
grid
hold off

% Exercise 5

clear
close all

m = 2
R = 1
k = 3
M = 4
b = 5

A = @(x) (m .* R .* (x .^ 2)) ./ sqrt(((k - M .* (x .^ 2)) .^ 2) + ((b .* x) .^ 2))
w = (sqrt(2) * k) / sqrt(2 * M * k - b ^ 2)
A(w)

x = 0:0.1:10
f = (m .* R .* (x .^ 2)) ./ sqrt(((k - M .* (x .^ 2)) .^ 2) + ((b .* x) .^ 2))


figure(1)
hold on


%[fmax,indmax] = max(f) %largest value of vector f is its element number indmin
%xmax=x(indmax)
plot(x, f, 'linewidth',1.5)
ylabel('A(w)')
xlabel('w')

if(2 * M * k - b ^ 2 < 0)
  title(['m = ', num2str(m), ', R = ', num2str(R), ', k = ', num2str(k), ', M = ', num2str(M), ', b = ', num2str(b)])
else {
  plot(w, A(w), 'r.','markersize', 10)
  title(['m = ', num2str(m), ', R = ', num2str(R), ', k = ', num2str(k), ', M = ', num2str(M), ', b = ', num2str(b), ' : max A(', num2str(w), ')', ' = ', num2str(A(w))])
  }
end
grid
hold off

% POSITION VELOCITY ACCELERATION:

% Exercise 1

clear
close all
# GIVEN
x = 3
y = 4
v = 2


z = sqrt(x ^ 2 + y ^ 2)
%position
s = @(t) sqrt((z + 2 * t) .^ 2 - y .^ 2)
%velocity
v = @(t) (2 .* (2 .* t + z)) ./ sqrt((2 .* t + z) .^2 - y .^2)
%acceleration
a = @(t) -(4 .* y .^ 2) ./ ((2 .* t + z) .^ 2 - y .^ 2) .^ (3 ./ 2)
%values at time t0
t0=2
s(t0)
v(t0)
a(t0)


%graphs when t=0...T
T=2
t=0:T/100:T;


subplot(311)
plot(t,s(t),'linewidth',1)
grid
xlim([0,T])
title('position s(t)')
subplot(312)
plot(t,v(t),'linewidth',1)
grid
xlim([0,T])
title('velocity v(t)')
subplot(313)
plot(t,a(t),'linewidth',1)
grid
xlim([0,T])
title('acceleration a(t)')
xlabel('time t')

% Exercise 2

h = 5 %height OA
x = pi / 6 % angle


w = 2 * pi %angular velocity of P
a = -5 * pi %angular acceleration of P


c = tan(x)* h % position


u = w * (c^2+h^2)/ h % velocity


v = (((a * ((c ^ 2 + h ^ 2) ^ 2) / h) + 2 * c * (u ^ 2)) / (c ^ 2 + h ^ 2)) % acceleration


%compare

% Exercise 3

clear
close all
format long
%GIVEN

w = 2 .* pi
T = 2 .* pi ./ w
dt = T ./ 100
t = 0:dt:T

E_D = 400
O_A = sqrt(100 .^ 2 + 200 .^ 2)
O_D = 600
A_D = sqrt(100 .^ 2 + 400 .^ 2)
A_B = sqrt(300 .^ 2 + 50 .^ 2)
C_x = 200
C_y = 100
B_C = 50

% Solve
alfa = w .* t
% Coordinates of B
B_x = (C_x + sin(alfa) .* B_C)
B_y = (C_y + cos(alfa) .* B_C)

% OB and B
O_B = sqrt(B_x .^ 2 + B_y .^ 2)
B = acos(B_x ./ O_B)

% AOB and AOD
A_O_B = acos(((A_B .* A_B) - (O_A .* O_A) - (O_B .* O_B)) ./ (- 2 .* O_A .* O_B))
A_O_D = (acos(((A_D .* A_D) - (O_A .* O_A) - (O_D .* O_D)) ./ (- 2 .* O_A .* O_D))) % h

% polar angle
polar_angle = (A_O_B - A_O_D + B)

u = sin(polar_angle)
u_c =  cos(polar_angle)

% Final result
D_x = O_D .* u_c
E_x = (D_x + sqrt(E_D .^ 2 - ((O_D .* u - 400) .^ 2))) / 1000

% graph
figure(1)
n = length(w .* t);
vk = (E_x(2:n) - E_x(1:n-1)) ./ dt %n-1 numbers
ak=(vk(2:n-1)-vk(1:n-2))/dt; %n-2 numbers

subplot(3, 1, 1)
plot(t, E_x, 'linewidth', 1.5)
title('E:n paikka s = Ex (m)')
grid
subplot(3, 1, 2)
plot(t(1:n-1),vk, 'linewidth', 1.5)
grid
title('E:n nopeus v (m / s)')
subplot(3,1,3)
plot(t(1:n-2),ak,'linewidth',1.5)
grid
title('E:n kiihtyvyys a (m / s2)')
xlabel('aika t (s)')
hold off

% ACCELERATION VELOCITY POSITION

% Exercise 1

clear
close all

%GIVEN
T1 = 1
T2 = 2
T3 = 3
amax = 2
T = 4 * T1 + 2 * T2 + T3

t1 = T1
t2 = T1 + T2
t3 = T1 + T2 + T1
t4 = T1 + T2 + T1 + T3
t5 = T1 + T2 + T1 + T3 + T1
t6 = T1 + T2 + T1 + T3 + T1 + T2
t7 = T1 + T2 + T1 + T3 + T1 + T2 + T1


dt = T / 1000
t = 0:dt:T;


a1 = amax * t
a2 = amax
a3 = -amax * t + 8
a4 = 0
a5 = -amax * t + 14
a6 = -amax
a7 = amax * t - 22

N=length(t)
vk=zeros(1,N);
sk=zeros(1,N);
vk(1)=0;
sk(1)=0;

a = a1.*(t<=t1) + a2.*((t>t1)&(t<t2)) + a3.*((t>t2)&(t<t3)) + a4.*((t>t3)&(t<t4)) + a5.*((t>t4)&(t<t5)) + a6.*((t>t5)&(t<t6)) + a7.*((t>t6)&(t<t7));

for n=1:N-1
    vk(n+1)=vk(n)+a(n)*dt;
    sk(n+1)=sk(n)+vk(n)*dt+1/2*a(n)*dt^2;
end

figure(4)
subplot(3,1,1)
plot(t,a,'linewidth',1.5)
xlim([0, T])
grid
title('kiihtyvyys a(t)')
subplot(3,1,2)
plot(t,vk,'linewidth',1.5)
xlim([0, T])
title('nopeus v(t)')
grid
subplot(3,1,3)
plot(t,sk,'linewidth',1.5)
xlim([0, T])
grid
title('paikka s(t)')
xlabel('aika t')

% Exercise 2

clear
close all

% GIVEN
m = 1
v_0 = 5
alpha = 65 / 180 * pi
b = 1
delta_t = 0.0001
g = 9.81
%SOLVE

%1) air resistance F(t)=0

dt = 0.01 %time step
s1x(1) = 0 %s(0)
s1y(1) = 0 %s(0)
v1x(1) = v_0 * cos(alpha) %v(0)
v1y(1) = v_0 * sin(alpha) %v(0)
k = 1
while s1y(k) >= 0
      F1x = 0 %air resistance
      F1y = 0
      a1x(k) = 1 / m * F1x;%acceleration
      a1y(k) = 1 / m * (F1y - m * g)

      v1x(k + 1) = v1x(k) + a1x(k) * dt;
      v1y(k + 1) = v1y(k) + a1y(k) * dt;

      s1x(k + 1) = s1x(k) + v1x(k) * dt + 1/2 * a1x(k) * dt ^ 2;
      s1y(k + 1) = s1y(k) + v1y(k) * dt + 1/2 * a1y(k) * dt ^ 2;
      k = k+1;
end

t1=0:dt:(k-1)*dt;
T1=t1(end)
s1=s1x(end)

%2)  F(t)=b*v

dt = 0.01 %time step
s2x(1) = 0 %s(0)
s2y(1) = 0 %s(0)
v2x(1) = v_0 * cos(alpha) %v(0)
v2y(1) = v_0 * sin(alpha) %v(0)
k = 1
while s2y(k) >= 0
      F2x = -b * v2x(k) %air resistance
      F2y = -b * v2y(k)
      a2x(k) = 1 / m * F2x;%acceleration
      a2y(k) = 1 / m * (F2y - m * g)

      v2x(k + 1) = v2x(k) + a2x(k) * dt;
      v2y(k + 1) = v2y(k) + a2y(k) * dt;

      s2x(k + 1) = s2x(k) + v2x(k) * dt + 1/2 * a2x(k) * dt ^ 2;
      s2y(k + 1) = s2y(k) + v2y(k) * dt + 1/2 * a2y(k) * dt ^ 2;
      k = k+1;
end

t2=0:dt:(k-1)*dt;
T2=t2(end)
s2=s2x(end)


%3) F(t)=b*v^2

dt = 0.01 %time step
s3x(1) = 0 %s(0)
s3y(1) = 0 %s(0)
v3x(1) = v_0 * cos(alpha) %v(0)
v3y(1) = v_0 * sin(alpha) %v(0)
k = 1
while s3y(k) >= 0
      F3x = -b * v3x(k) * v3x(k) %air resistance
      F3y = -b * v3y(k) * v3y(k)
      a3x(k) = 1 / m * F3x;%acceleration
      a3y(k) = 1 / m * (F3y - m * g)

      v3x(k + 1) = v3x(k) + a3x(k) * dt;
      v3y(k + 1) = v3y(k) + a3y(k) * dt;

      s3x(k + 1) = s3x(k) + v3x(k) * dt + 1/2 * a3x(k) * dt ^ 2;
      s3y(k + 1) = s3y(k) + v3y(k) * dt + 1/2 * a3y(k) * dt ^ 2;
      k = k + 1;
end

t3=0:dt:(k-1)*dt;
T3=t3(end)
s3=s3x(end)

figure(1)
hold on
plot(s1x,s1y,'r','linewidth',1.5)
plot(s2x,s2y,'green','linewidth',1.5)
plot(s3x,s3y,'blue','linewidth',1.5)

grid
title(['m = ',num2str(m),', v_0 = ',num2str(v_0),', a = ',num2str(alpha), ', b = ',num2str(b), ', alpha = ',num2str(65), ', delta_t = ',num2str(delta_t),
      'T_1 = ',num2str(T1), ', T_2 = ',num2str(T2), ', T_3 = ',num2str(T3), ', s_1 = ',num2str(s1), ', s_2 = ',num2str(s2), ', s_3 = ',num2str(s3),
])


% Exercise 3

clear
close all

%GIVEN
L1 = 0.5
L2 = 1

%SOLVE

V = 1 %speed of front wheel
T = 10 % End time
dt = T / 1000 % time step
t = 0:dt:T
N = length(t)

%trajectory of joint wheel
x1 = zeros(1,N) % zero from 1 to N
y1 = zeros(1,N)

%1) front wheel moves to the right on x-axis with speed V

X = V*t; %X(t)
Y = zeros(1,N); %Y(t)=0
X(1) = 0
Y(1) = 0
VX = V*ones(1,N); %X'(t)=V
VY = zeros(1,N); %Y'(t)=0

x1(1) = 0
y1(1) = L1

%2) front wheel moves on a circle with speed V
%r=2
%X=r*cos(V*t/r);
%Y=r*sin(V*t/r);
%VX=-V*sin(V*t/r);%X'(t)
%VY=V*cos(V*t/r);%Y'(t)
% %starting point of the rear wheel
%x1(1)=r;
%y1(1)=-L1;

for k = 1:N-1
   vx1(k) = ((X(k)-x1(k))*VX(k)+(Y(k)-y1(k))*VY(k))/L1^2*(X(k)-x1(k)); %x'(t)
   vy1(k) = ((X(k)-x1(k))*VX(k)+(Y(k)-y1(k))*VY(k))/L1^2*(Y(k)-y1(k)); %y'(t)
   dx1 = vx1(k) * dt;
   dy1 = vy1(k) * dt;
   x1(k + 1) = x1(k) + dx1;
   y1(k + 1) = y1(k) + dy1;
end

%2 takarengas
%2) joint wheel moves to the right on x-axis with speed V

X1 = V * t
Y1 = zeros(1,N); %Y(t)=0

X1(1) = 0
Y1(1) = L1

VX1 = V*ones(1,N); %X'(t)=V
VY1 = zeros(1,N); %Y'(t)=0

%trajectory of rear wheel
x2 = zeros(1,N) % zero from 1 to N
y2 = zeros(1,N)

x2(1) = -1
y2(1) = L1

%2) front wheel moves on a circle with speed V
%r = 2
%X1 = r * cos(V * t / r);
%Y1 = r * sin(V * t / r);
%VX1 = -V * sin(V * t / r);%X'(t)
%VY1 = V * cos(V * t / r);%Y'(t)
% %starting point of the rear wheel
%X1(1)=r;
%Y1(1)=-L1;

%x2(1) = 3
%y2(1) = -L1

for k = 1:N-1
   vx2(k) = ((X1(k)-x2(k))*VX1(k)+(Y1(k)-y2(k))*VY1(k))/L2^2*(X1(k)-x2(k)); %x'(t)
   vy2(k) = ((X1(k)-x2(k))*VX1(k)+(Y1(k)-y2(k))*VY1(k))/L2^2*(Y1(k)-y2(k)); %y'(t)
   dx2 = vx2(k) * dt;
   dy2 = vy2(k) * dt;
   x2(k + 1) = x2(k) + dx2;
   y2(k + 1) = y2(k) + dy2;
end



%% animation

figure(1)
for k = 1: N
plot(x2,y2,'b','linewidth',1.5)
%k=200
hold on
plot(X,Y,'r','linewidth',1.5)
plot([x1(k),X(k)],[y1(k),Y(k)],'k-','linewidth',2)
plot(X(k),Y(k),'r.','markersize',20)
plot(x1(k),y1(k),'black.','markersize',20)

plot([x2(k),x1(k)],[y2(k),y1(k)],'k-','linewidth',2)

plot(x2(k),y2(k),'blue.','markersize',20)

hold off
grid
axis equal
pause(0.01)
end
title(['L1 = ',num2str(L1),', L2 = ',num2str(L2)])

% Exercise 4

clear
close all

% GIVEN


V = 1 % speed of target
v = 1.1 % speed of the pursuer

T = 10 %end time
dt = T / 1000 %time step
t = 0:dt:T;
N = length(t);

% SOLVE

% 1) target moves to the right on x-axis with speed V
X = V * t; %X(t) = Vt
Y = zeros(1, N); %Y(t) = 0

VX = V * ones(1, N); %X'(t)=V
VY = zeros(1, N); %Y'(t)=0

%trajectory of the pursuer
x = zeros(1,N) % zero from 1 to N
y = zeros(1,N)

x(1) = 0
y(1) = 1


for k=1:N-1
   vx(k) = v * (X(k) - x(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
   vy(k) = v * (Y(k) - y(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
   dx = vx(k) * dt;
   dy = vy(k) * dt;
   x(k + 1) = x(k) + dx;
   y(k + 1) = y(k) + dy;
end

% 2) target moves on a circle

% GIVEN

%V = 1
%v = 1.1

%r = 1

%trajectory of the target
%X=r*cos(V*t/r);
%Y=r*sin(V*t/r);
%X(1) = r
%Y(1) = 0
%VX=-V*sin(V*t/r);%X'(t)
%VY=V*cos(V*t/r);%Y'(t)

%trajectory of the pursuer

%x = zeros(1,N) % zero from 1 to N
%y = zeros(1,N)

%x(1) = 0
%y(1) = 0

%for k=1:N-1
%   vx(k) = v * (X(k) - x(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
%   vy(k) = v * (Y(k) - y(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
%   dx = vx(k) * dt;
%   dy = vy(k) * dt;
%   x(k + 1) = x(k) + dx;
%   y(k + 1) = y(k) + dy;
%end


figure(1)
for k = 1:N
%k=200
plot(x,y,'b','linewidth',1.5)
hold on
plot(X,Y,'r','linewidth',1.5)
plot([x(k),X(k)],[y(k),Y(k)],'k-','linewidth',2)
plot(X(k),Y(k),'r.','markersize',20)
plot(x(k),y(k),'b.','markersize',20)
hold off
grid
axis equal
pause(0.01)
end

% SPLINE AND BEZIER CURVES:
% Exercise 1

clear
close all
x=[1, 4, 3, 1, 0]
y=[1, 2, 4, 5, 3]
n=length(x)

%tapaus=case in finnish
%tapaus=1 %1=natural, 2=clamped, 3=periodic

%tapaus=3 clamped spline
%k1=1 %Y'(x1)
%kn=2 %Y'(xn)

tapaus = 3
% then requirement y(1)=y(n)


%coefficients from a function spline_curves_coefficients.m
%which has to be on the current folder
%see spline_curves_coefficients.pdf
n=length(x)

%tapaus=3 %1=natural, 2=clamped, 3=periodic
%system of equations for  zk=Y''(xk)

A=zeros(n,n); %coefficient matrix
B=zeros(n,1); %right hand side
h=x(2:end)-x(1:end-1); %lenghts of the intervals x_k...x_(k+1)

for k=2:n-1
  alfak=h(k-1)/6;
  betak=(h(k-1)+h(k))/3;
  gammak=h(k)/6;
  deltak=(y(k+1)-y(k))/h(k)-(y(k)-y(k-1))/h(k-1);
  A(k,:)=[zeros(1,k-2),alfak,betak,gammak,zeros(1,n-(k-2)-3)];
  B(k)=deltak;
end


if tapaus==1
%1) natural spline

A(1,:)=[1,zeros(1,n-1)];
A(n,:)=[zeros(1,n-1),1];
B(1)=0;
B(n)=0;

elseif tapaus==2
%2) clamped
%k1=0 %Y'(x1)
%kn=-2 %Y'(xn)

alfa1=h(1)/3
beta1=h(1)/6
delta1=(y(2)-y(1))/h(1)-k1
alfan=h(n-1)/6
betan=h(n-1)/3
deltan=kn-(y(n)-y(n-1))/h(n-1)

A(1,:)=[alfa1,beta1,zeros(1,n-2)];
A(n,:)=[zeros(1,n-2),alfan,betan];
B(1)=delta1;
B(n)=deltan;

else
%3) periodic

alfa1=h(1)/3
beta1=h(1)/6
mu1=h(n-1)/6
la1=h(n-1)/3
delta1=(y(2)-y(1))/h(1)-(y(n)-y(n-1))/h(n-1)


A(1,:)=[alfa1,beta1,zeros(1,n-4),mu1,la1];
A(n,:)=[1,zeros(1,n-2),-1];
B(1)=delta1;
B(n)=0;
end

%solution i.e z1,z2,...,zn

z=A^-1*B

%coefficients for curves Yk=ak(x-xk)^3+bk(x-xk)^2+ck(x-xk)+dk, k=1,2,...,n-1,

for k=1:n-1
a(k)=(z(k+1)-z(k))/(6*h(k));
b(k)=z(k)/2;
c(k)=-h(k)/3*z(k)-h(k)/6*z(k+1)+(y(k+1)-y(k))/h(k);
d(k)=y(k);
end

abcd=[a;b;c;d]'

%Y_k=a_k*x^3+b_k*x^2+c_k*x+d_k
%coefficients a_k,b_k,c_k,d_k
%as rows in matrix abcd
a=abcd(:,1);
b=abcd(:,2);
c=abcd(:,3);
d=abcd(:,4);
%% graph
figure(1)
plot(x,y,'r.','markersize',18)
hold on
for k=1:n-1
  dx=(x(k+1)-x(k))/100;
  xx=x(k):dx:x(k+1);
  yy=a(k)*(xx-x(k)).^3+b(k)*(xx-x(k)).^2+c(k)*(xx-x(k))+d(k);
  plot(xx,yy,'linewidth',1.5)
end
hold off
grid
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
if tapaus==1
  title('natural spline','fontsize',10)
elseif tapaus==2
  title('clamped spline','fontsize',10)
else
  title('periodic spline','fontsize',10)
end

% INTEGRAL
% EXERCISE 1

clear
close all

format long

a = 0.7
h = 4

area = (2 * h * sqrt(h / a)) / 3

%Coordinates of the centroid
xcentroid = 1 ./ area .* -((h - 8) * h) / (4 * a)
%y
ycentroid = 1 ./ area .* (2 * h ^ 2 * sqrt(h / a)) / 5

% graph
hold on
plot([0, 0], [0, h], 'black','linewidth', 1.7)
plot([0, sqrt(h ./ a)], [h, h], 'b','linewidth', 1.7)


plot([-2, h], [0, 0], 'black','linewidth', 0.5)
plot([0, 0], [-1, 5], 'black','linewidth', 0.5)


x = 0:0.1:sqrt(h ./ a) + 0.1
f = @(x) (a .* (x .^ 2))
plot(x, f(x), 'r','linewidth', 1.4)
xlim([-2, 4])
ylim([-1, 5])
xlabel('x')
grid
title(['a = ', num2str(a) ,', h = ',num2str(h), ': A = ', num2str(area), ', xp = ', num2str(xcentroid), ', yp = ', num2str(ycentroid)])
axis equal
plot(xcentroid, ycentroid, 'black.' ,'markersize', 17)
hold off

% Exercise 2

clear
close all

% GIVEN

A = 5
T = 4

T = T / 2
% Cal
% u(t) = A / T * x

% Integrate THE MEAN
% The mean
uavg = 1 / T * A .* T ./ 2


% Integrate THE ROOT MEAN SQUARE
% THE ROOT MEAN SQUARE
urms = sqrt(1 / (T) * (A .^ 2) .* T / 3)

% Exercise 3

clear
close all
format long

% v1 = 500, p1 = 30
v_1 = 500, p_1 = 30
% v2 = 50, p2 = 475.4679577
v_2 = 50, p_2 = 475.4679577
% v3 = 50, p3 = 600
v_3 = 50, p_3 = 600
% v4 = 100, p4 = 600
v_4 = 100, p_4 = 600
% v5 = 200, p5 = 300
v_5 = 200, p_5 = 300
% v6 = 500, p6 = 57.65397283
v_6 = 500, p_6 = 57.65397283
% k1 = 1.2 k2 = 1.8
k_1 = 1.2, k_2 = 1.8



% v2 - v4  | p_3 - p1 * (v1 / v2) ^ k1
sum1 = (v_4^k_1*((k_1-1)*p_3*v_2^k_1*v_4+(1-k_1)*p_3*v_2^(k_1+1)-p_1*v_1^k_1*v_2)+p_1*v_1^k_1*v_2^k_1*v_4)/((k_1-1)*v_2^k_1*v_4^k_1)
% v4 - v5 | p4 * (v4 / v5) - p1 * (v1 / v2) ^ k1
sum2 = ((k_1-1)*p_4*v_4^(k_1+1)*v_5^k_1*log(v_5)+((1-k_1)*p_4*v_4^(k_1+1)*log(v_4)-p_1*v_1^k_1*v_4)*v_5^k_1+p_1*v_1^k_1*v_4^k_1*v_5)/((k_1-1)*v_4^k_1*v_5^k_1)
% v5 - v1 |  p5 * (v5 / v6) ^ k2 - p1 * (v1 / v2) ^ k1

sum3 = -((k_1-1)*p_5*v_1*v_5^(k_2+k_1)+v_5^k_1*((1-k_1)*p_5*v_1^k_2*v_5+(1-k_2)*p_1*v_1^(k_2+1))+(k_2-1)*p_1*v_1^(k_2+k_1)*v_5)/((k_1-1)*(k_2-1)*v_1^k_2*v_5^k_1)

% Area of the drawing
area = sum1 + sum2 + sum3

% Graph

figure(1)
hold on
x = v_5:1:v_6
y = @(x) p_5 .* (v_5 ./ x) .^ (k_2)
plot(x, y(x), 'b', 'linewidth', 1.5)

x_1 = v_2:1:v_6
y_1 = @(x_1) p_1 .* (v_1 ./ x_1) .^ (k_1)
plot(x_1, y_1(x_1), 'b', 'linewidth', 1.5)

x_2 = v_4:1:v_5
y_2 = @(x_2) p_4 .* (v_4 ./ x_2)
plot(x_2, y_2(x_2), 'b', 'linewidth', 1.5)


plot([v_1, v_6], [p_1, p_6], 'b', 'linewidth', 1.5)
plot([v_3, v_4], [p_3, p_4], 'b', 'linewidth', 1.5)
plot([v_2, v_3], [p_2, p_3], 'b', 'linewidth', 1.5)

plot(v_1, p_1, 'r.', 'markersize', 15)
plot(v_2, p_2, 'r.', 'markersize', 15)
plot(v_3, p_3, 'r.', 'markersize', 15)
plot(v_4, p_4, 'r.', 'markersize', 15)
plot(v_5, p_5, 'r.', 'markersize', 15)
plot(v_6, p_6, 'r.', 'markersize', 15)

xlabel('tilavuus V')
ylabel('paine p')
title(['v1 = ',num2str(v_1),', p1 = ',num2str(p_1), ', v2 = ', num2str(v_2), ', p3 = ', num2str(p_3), ', v4 = ', num2str(v_4), ', v5 = ', num2str(v_5), ', k1 = ', num2str(k_1), ', k2 = ', num2str(k_2), ', Ala = ', num2str(area)])
grid
hold off

% Exercise 4

clear
close all

% A(y) = sqrt(y) * 3/4 * (4 - y) dy

% GIVEN
h = 3
L = 4

% Area

% A = @(y) sqrt(y) .* h ./ L .* (L - y)
A = (4 * L ^ (3 / 2) * h) / 15

% Exercise 5

clear
close all
format long


% GIVEN
v_0 = 10
a = pi / 6
h = 5
g = 9.81

% Solve
% - g / (2 *(v0 * cos(a))^2) * (x^2) + tan(a) * x + h

first = - g / (2 * (v_0 * cos(a)) ^ 2)
second = tan(a)
delta = second ^ 2 - 4 * first * h
x1 = (- second + sqrt(delta)) / (2 * first)
x2 = (- second - sqrt(delta)) / (2 * first)


if(x1 > 0)
  k = x1
else
  k = x2
end

% The length s of the picture

s = -(2^(3/2)*cos(a)^4*v_0^4*asinh((sec(a)^2*sin(2*a)*v_0^2-2*sec(a)^2*g*k)/(2*v_0^2))+(2*cos(a)^2*tan(a)*v_0^2-2*g*k)*sqrt((cos(2*a)+1)*v_0^4-2*sin(2*a)*g*k*v_0^2+2*g^2*k^2)+(-2^(3/2)*cos(a)^4*asinh((sec(a)^2*sin(2*a))/2)-2*cos(a)^2*tan(a)*sqrt(cos(2*a)+1))*v_0^4)/(2^(5/2)*cos(a)^2*g*v_0^2)

% Graph
figure(1)
hold on
x = 0:0.1:k
f = @(x) -g ./ (2 .* (v_0 .* cos(a)) .^ 2) .* (x .^ 2) + tan(a) .* x + h
plot(x, f(x), 'linewidth', 1.5)
xlim([-1, k])
ylim([-2, v_0])
plot([0, k], [0, 0], 'black' ,'linewidth', 1)
plot([0, 0], [0, h + 3], 'black' ,'linewidth', 1)
title(['v_0 = ',num2str(v_0),', a = ',num2str(30), ', h = ', num2str(h), ', s = ', num2str(s)])

xlabel('x')
ylabel('y')
grid
hold off

% Exercise 6

clear
close all
format long


% t = 0 . . . 2π

%Given
R = 5
L = 3


% Derivative of x(t)
x = @(t) -3*L*sin(3*t)-R*sin(t)
% Derivative of y(t)
y = @(t) 3*L*cos(3*t)+R*cos(t)



N = 1000%number of time intervals
t1 = 0 %start time
t2 = 2*pi %end time
dt = (t2-t1)/N %time step

f = @(t) sqrt((x(t) .^ 2) + (y(t) .^ 2))

tk = (t1+dt/2):dt:t2; %center points of time intervals
vk = f(tk); %velocities
areas = vk*dt; %areas
area = sum(areas)

% graph

figure(1)
hold on
axis equal
t = 0:0.01:(2*pi)
u = @(t) R * cos(t) + L * cos(3 * t)
v = @(t) R * sin(t) + L * sin(3 * t)
plot(u(t),v(t),'blue' ,'linewidth', 1.5)
xlabel('x')
ylabel('y')
title(['R = ', num2str(R), ', L = ', num2str(L), ', s = ', num2str(area)])

grid
hold off

% Exercise 7

clear
close all

%Given
x0 = 0
y0 = 0

R = 2
a = 0.698

zero_l = 1.047

% Solve
zero_0 = 0
zero_l1 = (zero_l - a) / 2

k_0 = 0
k_l = 0

k_l1 = 1 / R

% x1, y1
% The length L of the curve is determined by
% L = 2(θ(L) − θ(0))  / (κ(L) + κ(0))


L1 = (2 * (zero_l1 - zero_0)) / ((k_l1) + (k_0))

ds=L1/100 %step
s=0:ds:L1;

N=length(s)
x=zeros(1,N);
y=zeros(1,N);

ths_1 = zero_0 + k_0 * s + (k_l1 - k_0)/(2*L1) * s .^ 2; %direction angle

x(1)=x0;
y(1)=y0;

for n=1:N-1;
   x(n+1)=x(n)+ds*cos(ths_1(n));
   y(n+1)=y(n)+ds*sin(ths_1(n));
end

figure(1)
hold on

plot(x,y,'r','linewidth',1.5)
plot(x(1),y(1),'r.','markersize',20)
q = x(end)
w = y(end)

xlim([-0.5, 2.5])
ylim([-1, 2])
grid

zero_l2 = zero_l1 + a
k_l2 = 1 / R

% line 2
L2 = (2 * (zero_l2 - zero_l1)) / ((k_l2) + (k_l1))

ds_2=L2/100 %step
s_2= 0:ds_2:L2;

N_2=length(s_2)
x_2=zeros(1,N_2);
y_2=zeros(1,N_2);

ths_2 = zero_l1 + k_l1 * s_2 + (k_l2 - k_l1) / (2 * L2) * s_2 .^ 2; %direction angle

x_2(1) = q;
y_2(1) = w;


for n_2=1:N_2-1;
   x_2(n_2+1)=x_2(n_2)+ds_2*cos(ths_2(n_2));
   y_2(n_2+1)=y_2(n_2)+ds_2*sin(ths_2(n_2));
end

e = x_2(end)
r = y_2(end)

plot(x_2,y_2,'black','linewidth',1.5)

% Line 3
zero_l3 = zero_l
k_l3 = 1 / R

L3 = (2 * (zero_l3 - zero_l2)) / ((k_l3) + (k_0))


ds_3 = L3/100 %step
s_3 = 0:ds_3:L3;

N_3=length(s_3)
x_3=zeros(1,N_3);
y_3=zeros(1,N_3);

ths_3 = zero_l2 + k_0 * s_3 + (k_l3 - k_0) / (2 * L3) * s_3 .^ 2; %direction angle

x_3(1) = e;
y_3(1) = r;

for n_3=1:N_3-1;
   x_3(n_3+1)=x_3(n_3)+ds_3*cos(ths_3(n_3));
   y_3(n_3+1)=y_3(n_3)+ds_3*sin(ths_3(n_3));
end

plot(x_3,y_3,'blue','linewidth',1.5)
plot(x_3(end),y_3(end),'blue.','markersize',20)

title([
     '\theta(L) = ',num2str(zero_l),...
     ', R = ',num2str(R),...
     ', a = ',num2str(a),
     '\theta(L1) = ',num2str(zero_l1),...
     ', L1 = ',num2str(L1),...
     ', L2 = ',num2str(L2),...
     ', L3 = ',num2str(L3),...
           ])
hold off

figure(2)
hold on

ks=k_0+(k_l1 - k_0)/L1*s; %curvature
plot(s, ks, 'r','linewidth',1.5)

s_22 = 0:ds_2:L2;


ks_2 = k_l1 + (k_l2 - k_l1) / L2 * s_22; %curvature

plot(s_22 + L1, ks_2, 'black','linewidth',1.5)


s_33 = 0:ds_3:L3;

ks_3 = k_l3 + (k_0 - k_l3) / L3 * s_33; %curvature
plot(s_33 + L1 + L2, ks_3, 'blue', 'linewidth', 1.5)

title([
     'kaarevuus k(s)'
           ])

grid
hold off

figure(3)
hold on
ths=zero_0 + k_0 * s + (k_l1 - k_0)/(2 * L1) * s.^2; %direction angle
plot(s, ths, 'r', 'linewidth', 1.5)

ths22 =  k_l1 * s_22 + ((k_l2 - k_l1) / ((2 * L2)) * ((s_22) .^ 2)); %direction angle
plot(s_2 + L1, ths_2, 'black', 'linewidth', 1.5)


ths33 = zero_l2 + k_l3 * s_33 + ((k_0 - k_l3) / (2 * L3)) * s_33.^2; %direction angle
plot(s_33 + L1 + L2, ths33, 'blue', 'linewidth', 1.5)
title([
     'suuntakulma ', '\theta(s)'
           ])
xlabel('käyrän pituus s')
grid
hold off

% TAYLOR POLYNOMIAL

% Exercise 1

clear
close all
#GIVEN

n = 5

#SOLVE

Tn = 0
x = -1:0.01:1
for k = 1:n
   Tn = Tn + (-1)  .^ (k - 1) ./ k .* x .^ (k)
end

#graph
hold on
f = @(x) log(x + 1)
plot(x, f(x), 'blue' ,'linewidth', 1.5)
plot(x,Tn,'r','linewidth',1.2)
grid
ylim([-2, 2])

title(['n = ', num2str(n)])
xlabel('x')

hold off

