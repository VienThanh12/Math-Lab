%% ex 1, p.3
clear
close all
%x'(t)+6/1000*x(t)=6, x(0)=0
%first degree equation
%a*y'(t)+b*y(t)=A,y(0)=y0
%solution_formulas.pdf,  1 b i), p.2
a=1
b=6/1000
A=6
y0=0
la=-b/a
%solution
x=@(t) A/b+(y0-A/b)*exp(la*t)

tau=-1/la %time constant
tmax=10*tau
t=0:tau/100:tmax;

plot(t,x(t),'linewidth',1.5)
grid
xlabel('time (min)')
title('mass of salt (kg)')
xlim([0,tmax])
%% ex 2, RC-circuit, p.5
clear
close all
R=1
C=2
U=4
%R*q'(t)+1/C*q(t)=U, q(0)=0
%first degree equation
%a*y'(t)+b*y(t)=A,y(0)=y0
%solution_formulas.pdf, 1b i) p.2
a=R
b=1/C
A=U
y0=0
la=-b/a
%solution
q=@(t) A/b+(y0-A/b)*exp(la*t)
%i(t)=q'(t)
i=@(t) (y0-A/b)*exp(la*t)*la

tau=-1/la %time constant
tmax=5*tau
t=0:tau/100:tmax;
plot(t,q(t),'linewidth',1.5)
hold on
plot(t,i(t),'linewidth',1.5)
hold off
grid
xlabel('time t')
xlim([0,tmax])
legend({'q(t)','i(t)'},'fontsize',12)

%% ex 2, RC-circuit, p.7
R=1
C=2
U=1
w=2*pi
%R*q'(t)+1/C*q(t)=U*sin(w*t), q(0)=0
%first degree equation
%a*y'(t)+b*y(t)=A*sin(w*t+theta)+B,y(0)=y0
%solution_formulas.pdf, 1 b ii), p.2
a=R
b=1/C
A=U
B=0
y0=0 %q(0)
la=-b/a
K=1/sqrt(b^2+(a*w)^2)
phi=atan2(-a*w,b)
D=y0-K*A*sin(phi)-B/b

q=@(t) K*A*sin(w*t+phi)+B/b+D*exp(la*t)

dq=@(t) K*A*cos(w*t+phi)*w+D*exp(la*t)*la %q'(t)=i(t)

u=@(t) A*sin(w*t)+B

tau=-1/la %time constant
tmax=5*tau

t=0:tau/100:tmax;

figure(1)
plot(t,q(t),'linewidth',1.5)
hold on
plot(t,u(t),'linewidth',1.5)
plot(t,K*A*sin(w*t+phi),'linewidth',1.5)
hold off
grid
xlabel('time t')
legend({'q(t)','u(t)','KUsin(\omega t+\phi)'},'fontsize',12)
title(['R = ',num2str(R),', C = ',num2str(C),', U = ',num2str(U),', \omega = ',num2str(w),', K = ',num2str(K),', \phi = ',num2str(phi)])


%check
t=1.432
a*dq(t)+b*q(t)
u(t)
q(0)
y0

%% amplification and phase shift
w=0:0.01:10;
K=1./sqrt(b^2+(a*w).^2)
phi=atan2(-a*w,b)

figure(2)
subplot(2,1,1)
plot(w,K,'linewidth',1.5)
grid
title('amplification K')
subplot(2,1,2)
plot(w,phi,'linewidth',1.5)
grid
title('phase shift \phi')
xlabel('angular frequency \omega')
%% Newtons law of cooling, p.13
clear
k=2
T0=10
T1=5
T2=20
w=10

%1*T'(t)+k*T(t)=k*T1*sin(w*t)+k*T2, T(0)=T0
%first degree equation
%a*y'(t)+b*y(t)=A*sin(w*t+theta)+B,y(0)=y0
%solution_formulas.pdf, 1 b ii) p.2

a=1
b=k
A=k*T1
B=k*T2
y0=T0
la=-b/a
K=1/sqrt(b^2+(a*w)^2)
phi=atan2(-a*w,b)
C=y0-K*A*sin(phi)-B/b


T3=K*A

T=@(t) K*A*sin(w*t+phi)+B/b+C*exp(la*t)

Ty=@(t) T1*sin(w*t)+T2

tau=-1/la %time constant
tmax=6*tau
t=0:tau/1000:tmax;

figure(1)
plot(t,T(t),'linewidth',1.5)
hold on
plot(t,Ty(t),'linewidth',1.5)
plot(t,T3*sin(w*t+phi)+T2,'linewidth',1.5)
hold off
grid
xlabel('time t')
legend({'T(t)','T_y(t)','T_3sin(\omega t+\phi)+T_2'},'fontsize',12)
title(['T_3 = ',num2str(T3),', \phi = ',num2str(phi)])
xlim([0,tmax])

%% amplitude T3 and phase shift phi
wmax=50
w=0:wmax/100:wmax;
K=1./sqrt(b^2+(a*w).^2);
T3=K*A;
phi=atan2(-a*w,b);

figure(2)
subplot(2,1,1)
plot(w,T3,'linewidth',1.5)
grid
title('amplitude T_3')
subplot(2,1,2)
plot(w,phi,'linewidth',1.5)
grid
title('phase shift \phi')
xlabel('angular frequency \omega')



%% damped oscillation, p.19
clear
close all
m=2
b=1
k=10
%b=2*sqrt(m*k)
x0=1
v0=2
%m*x''(t)+b*x'(t)+k*x(t)=0, x(0)=x0, x'(0)=v0
%second degree equation
%a*y''(t)+b*y'(t)+c*y(t)=0, y(0)=y0, y'(0)=y1
%solution_formulas.pdf 2 a), p.3
a=m
%b=b
c=k
y0=x0
y1=v0


%function second_degree_homogeneous.m
[x,dx,d2x,tau]=second_degree_homogeneous(a,b,c,y0,y1)

tmax=5*tau
t=0:tau/100:tmax;

figure(1)
subplot(3,1,1)
plot(t,x(t),'linewidth',1.5)
grid
xlabel('time t')
title('position x(t)')
xlim([0,tmax])


subplot(3,1,2)
plot(t,dx(t),'linewidth',1.5)
title('velocity v(t)')
grid
xlim([0,tmax])

subplot(3,1,3)
plot(t,d2x(t),'linewidth',1.5)
title('acceleration a(t)')
grid
xlabel('time t')
xlim([0,tmax])

%check
t=1.417
a*d2x(t)+b*dx(t)+c*x(t)
x(0)
y0
dx(0)
y1

%% forced oscillation, p.23-
clear
m=1
b=1
k=10
x0=0
v0=0
A=3
w=10*pi

%m*x''(t)+b*x'(t)+k*x(t)=A*sin(w*t), x(0)=x0, x'(0)=v0
%second degree equation
%a*y''(t)+b*y'(t)+c*y(t)=A*sin(w*t+theta)+B, y(0)=y0, x'(0)=y1
%solution_formulas.pdf 2 b ii), p.6
a=m
%b=b
c=k
y0=x0
y1=v0
B=0



K=1/sqrt((b*w)^2+(c-a*w^2)^2)
phi=atan2(-b*w,c-a*w^2)

xp=@(t) K*A*sin(w*t+phi)+B/c


y0h=y0-K*A*sin(phi)-B/c
y1h=y1-K*A*w*cos(phi)

[xh,dxh,d2xh,tau]=second_degree_homogeneous(a,b,c,y0h,y1h)

x = @(t) xp(t)+xh(t)

dx=@(t) K*A*cos(w*t+phi)*w +dxh(t)
d2x=@(t) -K*A*sin(w*t+phi)*w^2 +d2xh(t)

F=@(t) A*sin(w*t)

tmax=6*tau
t=0:tmax/1000:tmax;

figure(1)
subplot(2,1,1)
plot(t,A*sin(w*t),'linewidth',1.5)
grid
legend({'F(t)'},'fontsize',11)
title(['m = ',num2str(m),', b = ',num2str(b),', k = ',num2str(k),...
      ', A = ',num2str(A),', \omega = ',num2str(w),', x_0 = ',num2str(x0),', v_0 = ',num2str(v0)])
subplot(2,1,2)
plot(t,x(t),t,xp(t),'linewidth',1.5)
grid
xlabel('time t')
title(['K = ',num2str(K),', \phi = ',num2str(phi)])
legend({'x(t)','KAsin(\omega t +\phi)'},'fontsize',11)
%check
t=1.417
a*d2x(t)+b*dx(t)+c*x(t)
F(t)
x(0)
y0
dx(0)
y1

%% amplification K

if b<sqrt(2*k*m)
  w0=sqrt(k/m-b^2/(2*m^2))
  K0=1/(b*sqrt(k/m-b^2/(4*m^2)))
else
  w0=0
  K0=1/k
end

wmax=max([5*w0,5])
w=0:wmax/1000:wmax;


K=1./sqrt((b*w).^2+(c-a*w.^2).^2);
phi=atan2(-b*w,c-a*w.^2);

figure(2)
subplot(2,1,1)
plot(w,K,w0,K0,'r.','linewidth',1.5,'markersize',15)
grid
title(['amplification K: w_0 = ',num2str(w0),', K_0 = ',num2str(K0)])

subplot(2,1,2)
plot(w,phi,'linewidth',1.5)
grid
title('phase shift \phi')
xlabel('angular frequency t')

%% p.29, DC-motor
clear
close all
J=0.01
bm=0.1 %motor viscous friction constant
K=0.01
R=1
L=0.5
U=1%u(t)=U

%JL*w''(t)+(RJ+bmL)*w'(t)+(bR+K^2)*w(t)=K*u(t),w(0)=0,w'(0)=0
%second degree equation
%a*y''(t)+b*y'(t)+c*y(t)=A, y(0)=y0, x'(0)=y1
%solution_formulas.pdf, 2 b i), p.5

a=J*L
b=R*J+bm*L
c=bm*R+K^2
A=K*U
y0=0
y1=0

y0h=y0-A/c
y1h=y1
[wh,dwh,d2wh,tau]=second_degree_homogeneous(a,b,c,y0h,y1h)
w=@(t) A/c+wh(t)

tmax=10*tau
t=0:tmax/1000:tmax;

plot(t,w(t),'linewidth',1.5)
grid
xlabel('time t')
title('angular speed \omega(t)')

%w(t)-> A/c, as t increases
A/c


%% ex 7, p.33 RC-circuit, u(t)= square wave

clear
close all
R=1
C=0.5
q0=0
U=3%height of u(t)
T=2%period of u(t)

N=1000
dt=T/N;%time step
t=0:dt:T-dt;

%u(t), u=[U,U,...,U,0,0,...,0]
u=[U*ones(1,N/2),zeros(1,N/2)];
%ones(1,N) =[1,1,...,1] = N ones  (1xN-matrix)
%zeros(1,N)=[0,0,...,0] = N zeros   (1xN-matrix)

%copies for 5 period
u=[u,u,u,u,u];
t=[t,t+T,t+2*T,t+3*T,t+4*T];

%values of q(t) and i(t)=q'(t) numerically
N=length(t);
q=zeros(1,N);
dq=zeros(1,N-1);
q(1)=q0

for k=1:N-1
     dq(k)=(u(k)-1/C*q(k))/R; %q'(tk)=i(tk)
     q(k+1)=q(k)+dq(k)*dt;
end

subplot(2,1,1)
plot(t,u,'linewidth',1.5)
grid
ylim([-1 U+1])
title('voltage u(t)')

subplot(2,1,2)
plot(t,q,'linewidth',1.5)
hold on
plot(t(1:end-1),dq,'linewidth',1.5)
hold off
grid
legend({'q(t)','i(t)'},'fontsize',12)
xlabel('time t')
title('charge q(t) and current i(t)')
%% ex 8, p.35, mechanical oscillator, force F(t) = sawtooth
clear
close all
m=1
b=100
k=100
F0=100
T=2 %period
x0=0 %x(0)
v0=0 %x'(0)=v(0)


dt=0.001 %time step
%increasing part
t1=0:dt:T/2;
F1=F0/(T/2)*t1;
%decreasing part
t2=T/2+dt:dt:T-dt;
F2=-F0/(T/2)*(t2-T/2)+F0;
%one sawtooth
t=[t1,t2];
F=[F1,F2];

%5 sawtooth
t=[t,t+T,t+2*T,t+3*T,t+4*T];
F=[F,F,F,F,F];

figure(1)
subplot(2,1,1)
plot(t,F,'linewidth',1.5)
grid
title('force F(t)')

%% values of x(t),x'(t),x''(t) numerically
N=length(t)
x=zeros(1,N);
v=zeros(1,N);
a=zeros(1,N-1);
x(1)=x0;
v(1)=v0;


for n=1:N-1
   a(n)=1/m*(F(n)-b*v(n)-k*x(n)); %x''(t)=a(t)
   x(n+1)=x(n)+v(n)*dt+1/2*a(n)*dt^2;
   v(n+1)=v(n)+a(n)*dt;%x'(t)=v(t)
end


subplot(2,1,2)
plot(t,x,'linewidth',1.5)
grid
title('position x(t)')

%%
figure(2)
subplot(2,1,1)
plot(t,v,'linewidth',1.5)
grid
title('velocity  v(t)')

subplot(2,1,2)
plot(t(1:N-1),a,'linewidth',1.5)
grid
xlabel('time t')
title('acceleration a(t)')




%% ex 9, p.37, PID-controller
%differential equation of the system
%ay'(t)+by(t)=u(t),y(0)=0
%y'(t)=1/a*(u(t)-b*y(t))
clear
close all
a=2
b=5
Kp=10
Ki=5
Kd=0.1
dt=0.01 %time step
N=3000
t=0:dt:N*dt;
r=ones(1,N+1);%setpoint
y=zeros(1,N+1);
u=zeros(1,N);
e=zeros(1,N);
y(1)=0
integral=0
previous_error=0

for k=1:N
    error=r(k)-y(k);
    e(k)=error;
    integral=integral+error*dt;
    derivative=(error-previous_error)/dt;
    u(k)=Kp*error+Ki*integral+Kd*derivative;%output of PID
    dy=1/a*(u(k)-b*y(k));%y'(t)
    y(k+1)=y(k)+dy*dt;
    previous_error=error;
end


figure(1)
subplot(3,1,1)
plot(t(1:end-1),e,'linewidth',1.5)
grid
title('error e(t) = r(t)-y(t)')

subplot(3,1,2)
plot(t(1:end-1),u,'linewidth',1.5)
grid
title('u(t) = ouput of PID ')

subplot(3,1,3)
plot(t,y,'linewidth',1.5)
grid
xlabel('time t')
title('y(t) = output of the system')




