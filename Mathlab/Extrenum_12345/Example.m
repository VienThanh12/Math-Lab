%% p.3-4
clear
close all
x=0:0.01:3;
f=-x.^3+4*x.^2-3*x;
[fmin,indmin]=min(f)%smallest value of vector f is its element number indmin
xmin=x(indmin)
[fmax,indmax]=max(f) %largest value of vector f is its element number indmin
xmax=x(indmax)
%%
plot(x,f,'linewidth',1.5)
hold on
plot(xmin,fmin,'r.','markersize',20);
plot(xmax,fmax,'g.','markersize',20);
hold off
grid
xlabel('x')
ylabel('f(x)')

%% p.5-8
clear
close all
A=19
B=5

%minimum at x=x0
x0=sqrt(A/B)

%minimum value
L0=2*sqrt(A*B)

%graph
x=0:x0/100:3*x0;
L=A./x+B*x;

plot(x,L,'linewidth',1.5)
hold
plot(x0,L0,'r.','markersize',20)
hold off
grid
ylim([0,2*L0])
xlabel('x')
ylabel('L(x)')

%% p.9-10
clear
close all
H=3
h=1

%maximum at x=x0
x0=sqrt(H*h)
%maximum value
alfa0=atan(sqrt(H/h))-atan(sqrt(h/H))

%graph
x=0:x0/100:10*x0;
alfa=atan(H./x)-atan(h./x);

plot(x,alfa,'linewidth',1.5)
hold
plot(x0,alfa0,'r.','markersize',20)
hold off
grid
xlabel('x')
ylabel('\alpha(x)')

%% p.11-14
clear
close all
a=12
b=4
V=@(x) (a-2*x).*(b-2*x).*x
%maximum at
xmax=(a+b-sqrt((a-b)^2+a*b))/6
%maximum value
Vmax=V(xmax)
%Vmax=(-2*a^3 + a^2*(3*b + 2*sqrt(a^2 - a*b + b^2)+ a*b*(3*b - 2*sqrt(a^2 - a*b + b^2)) ...
%      + 2*b^2*(-b + sqrt(a^2 - a*b + b^2))))/54

%graph
x=0:b/100:b/2;
plot(x,V(x),'linewidth',1.5)
hold
plot(xmax,Vmax,'r.','markersize',20)
hold off
grid
xlabel('x')
ylabel('V(x)')

%% p.15-20

clear

k0=sqrt(29)/5%=1.077
k=1.02
x=0:0.01:10;
H=k*sqrt(x.^2+4^2)+1*(10-x);

figure(1)
if k<=k0
    x0=10 %minimum at x=x0
    H0=2*sqrt(29)*k %minimum value
    plot(x,H,'linewidth',1.5)
    hold
    plot(x0,H0,'r.','markersize',20)
    hold off
    grid
    xlabel('x')
    ylabel('H(x)')
    title(['k = ',num2str(k)])
else
    x0=4/sqrt(k^2-1) %minimum at x=x0
    H0=4*sqrt(k^2-1)+10 %minimum value
    plot(x,H,'linewidth',1.5)
    hold
    plot(x0,H0,'r.','markersize',20)
    hold off
    grid
    xlabel('x')
    ylabel('H(x)')
    title(['k = ',num2str(k)])
end

% geometry
figure(2)
plot([0,x0],[-4,0],'r','linewidth',2)
hold on
plot([x0,10],[0,0],'g','linewidth',2)
plot(x0,0,'b.','markersize',20)
plot([0,x0],[0,0],'k')
hold off
grid
axis equal
title(['x_0 = ',num2str(x0)])

%% p.21-26
clear
close all
k=10
m=10
b=10

D=b^2-2*k*m

if b^2-2*k*m>0 %maximum at w=0
    w=0:0.01:10;
    C=1./sqrt((k-m*w.^2).^2+(b*w).^2);
    plot(w,C,'linewidth',1.5)
    grid
    title('case 1')
else
    w0=sqrt(k/m-b^2/(2*m^2)); %maximum at w=w0
    C0=1/(b*sqrt(k/m-b^2/(4*m^2))); %maximum value
    w=0:w0/100:3*w0;
    C=1./sqrt((k-m*w.^2).^2+(b*w).^2);
    plot(w,C,'linewidth',1.5)
    hold
    plot(w0,C0,'r.','markersize',20)
    hold off
    grid
    xlim([0,3*w0])
    title('case 2')
end

%% p.27-32
clear
close all
L=7
h=2
H=3
v1=3
v2=5

x=0:L/1000:L;
AP=sqrt(x.^2+h^2);
PB=sqrt((L-x).^2+H^2);
T=AP/v1+PB/v2;

[Tmin,indmin]=min(T);
xmin=x(indmin);


%Snell's law
sin_alfa=xmin/sqrt(xmin^2+h^2)
sin_beta=(L-xmin)/sqrt((L-xmin)^2+H^2)
sin_alfa/v1
sin_beta/v2
alfa=asind(sin_alfa)
beta=asind(sin_beta)

figure(1)
plot(x,T,'linewidth',1.5)
hold
plot(xmin,Tmin,'r.','markersize',20)
hold off
grid
xlabel('x')
ylabel('T(x)  ','rotation',0)




%geometry
A=[0,-h]
P=[xmin,0]
B=[L,H]

figure(2)
plot([A(1),P(1),B(1)],[A(2),P(2),B(2)],'r.-','linewidth',2,'markersize',20)
hold on
plot([0,L],[0,0],'k')
plot([xmin,xmin],[-h,H],'k')
hold off
grid
axis equal
title(['v_1 = ',num2str(v1),', v_2 = ',num2str(v2),', xmin = ',num2str(xmin),...
       ', \alpha = ',num2str(alfa),', \beta = ',num2str(beta) ])



