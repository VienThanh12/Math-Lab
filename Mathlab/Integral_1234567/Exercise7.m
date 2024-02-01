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
