clear
close all

#Given
L = 0.5
R = 1
C = 2
U0 = 10
T = 2 % period

m = L
b = R
k = 1 / C
F0 = U0

x0 = 0 %x(0)
v0 = 0 %x'(0)=v(0)

dt = 0.001 %time step
%increasing part
t1 = 0:dt:T;  % 0:0.1:2.5
U1 = U0 / T * t1; % 10 / 5
%decreasing part
t2 = T + dt:dt:T - dt;
U2= -U0 / (T / 2) * (t2 - T / 2) + U0;
%one sawtooth
t = [t1,t2];
U = [U1,U2];

%5 sawtooth
t = [t, t + T, t + 2 * T, t + 3 * T, t + 4 * T];
U = [U, U, U, U, U];

figure(1)
subplot(2,1,1)
plot(t, U,'linewidth',1.5)
grid
ylabel('voltage u(t)')



N=length(t)
x=zeros(1,N);
v=zeros(1,N);
a=zeros(1,N-1);
x(1)=x0;
v(1)=v0;


for n=1:N-1
   a(n)=1/m*(U(n)-b*v(n)-k*x(n)); %x''(t)=a(t)
   x(n+1)=x(n)+v(n)*dt+1/2*a(n)*dt^2;
   v(n+1)=v(n)+a(n)*dt;%x'(t)=v(t)
end

figure(2)
subplot(3,1,1)
plot(t, x / C, 'red', 'linewidth', 1.5)
ylabel('voltage U_C')
grid
title(['L = ', num2str(L),', R = ', num2str(R), ', C = ', num2str(C), ', U = ', num2str(U0), ', T = ', num2str(T)])


%%
subplot(3,1,2)
plot(t,v * R, 'green','linewidth',1.5)
ylabel('voltage U_R')
grid

subplot(3,1,3)
plot(t(1:N-1),a * L, 'blue','linewidth',1.5)
ylabel('voltage U_L')
grid
xlabel('time t')




