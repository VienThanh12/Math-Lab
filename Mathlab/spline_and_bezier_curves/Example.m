clear
close all
x=[0,1,1.5,2.5,4,7]
y=[1,2,1.5,3,1,1]
n=length(x)

%tapaus=case in finnish
tapaus=1 %1=natural, 2=clamped, 3=periodic

%tapaus=2, clamped spline
k1=0 %Y'(x1)
kn=-2 %Y'(xn)

%if tapaus = 3, then requirement y(1)=y(n)


%coefficients from a function spline_curves_coefficients.m
%which has to be on the current folder
%see spline_curves_coefficients.pdf
abcd=spline_curves_coefficients(x,y,tapaus,k1,kn);
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
%% derivatives
figure(2)
subplot(2,1,1)
for k=1:n-1
  if k==2
      hold on
  end
  dx=(x(k+1)-x(k))/100;
  xx=x(k):dx:x(k+1);
  dyy=3*a(k)*(xx-x(k)).^2+2*b(k)*(xx-x(k))+c(k);
  plot(xx,dyy,'linewidth',1.5)
end
hold off

grid
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
title('derivative Y^\prime(x)')

subplot(2,1,2)

for k=1:n-1
  if k==2
      hold on
  end
  dx=(x(k+1)-x(k))/100;
  xx=x(k):dx:x(k+1);
  d2yy=6*a(k)*(xx-x(k))+2*b(k);
  plot(xx,d2yy,'linewidth',1.5)
end
hold off

grid
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
title('second derivative Y^{\prime\prime}(x)')


%% Cubic Bezier curve
clear
close all
%control points
P1=[0,0]
P2=[1,6]
P3=[4,5]
P4=[6,2]

x1=P1(1)
y1=P1(2)
x2=P2(1)
y2=P2(2)
x3=P3(1)
y3=P3(2)
x4=P4(1)
y4=P4(2)

x=@(t) (1-t).^3*x1+3*(1-t).^2.*t*x2+3*(1-t).*t.^2*x3+t.^3*x4;
y=@(t) (1-t).^3*y1+3*(1-t).^2.*t*y2+3*(1-t).*t.^2*y3+t.^3*y4;

t=0:0.01:1;

t0=0.6

figure(1)
plot([x1,x2,x3,x4],[y1,y2,y3,y4],'k.-','linewidth',1.5,'markersize',20)
hold
plot(x(t),y(t),'b','linewidth',1.5)
plot(x(t0),y(t0),'r.','markersize',20)
hold off
grid
axis equal
xlabel('x','fontsize',14)
ylabel('y  ','rotation',0,'fontsize',14)
title(['t_0 = ',num2str(t0)])
%% curvature
dx=-3*(1-t).^2*x1+...
   3*(1-4*t+3*t.^2)*x2+...
   3*(2*t-3*t.^2)*x3+...
   3*t.^2*x4;
d2x=6*(1-t)*x1+...
   3*(-4+6*t)*x2+...
   3*(2-6*t)*x3+...
   6*t*x4;
dy=-3*(1-t).^2*y1+...
   3*(1-4*t+3*t.^2)*y2+...
   3*(2*t-3*t.^2)*y3+...
   3*t.^2*y4;
d2y=6*(1-t)*y1+...
   3*(-4+6*t)*y2+...
   3*(2-6*t)*y3+...
   6*t*y4;

k=(dx.*d2y-d2x.*dy)./(dx.^2+dy.^2).^(3/2);

figure(2)
plot(t,k,'linewidth',1.5)
grid
xlabel('t')
ylabel('curvature \kappa(t)')

%% connecting Bezier curves
clear
close all
%control points
P1=[0,0]
P2=[1,3]
P3=[4,4]
P4=[5,2]

x1=P1(1)
y1=P1(2)
x2=P2(1)
y2=P2(2)
x3=P3(1)
y3=P3(2)
x4=P4(1)
y4=P4(2)


M=P3+(P3-P2) %=Q2+(Q2-Q3)

Q1=P4
Q2=Q1+(P4-P3)
%Q2=[2,-1]
Q3=Q2+(Q2-M)
%Q3=[4,-1]
Q4=[10,-2] %arbitrary

X1=Q1(1)
Y1=Q1(2)
X2=Q2(1)
Y2=Q2(2)
X3=Q3(1)
Y3=Q3(2)
X4=Q4(1)
Y4=Q4(2)


t=0:0.01:1;
x=(1-t).^3*x1+3*(1-t).^2.*t*x2+3*(1-t).*t.^2*x3+t.^3*x4;
y=(1-t).^3*y1+3*(1-t).^2.*t*y2+3*(1-t).*t.^2*y3+t.^3*y4;

X=(1-t).^3*X1+3*(1-t).^2.*t*X2+3*(1-t).*t.^2*X3+t.^3*X4;
Y=(1-t).^3*Y1+3*(1-t).^2.*t*Y2+3*(1-t).*t.^2*Y3+t.^3*Y4;

figure(2)

plot(x,y,'b','linewidth',1.5)
hold
plot(X,Y,'c','linewidth',1.5)
plot([x1,x2,x3,x4],[y1,y2,y3,y4],'k.-','linewidth',1.5,'markersize',20)
plot([X1,X2,X3,X4],[Y1,Y2,Y3,Y4],'.-','linewidth',1.5,'markersize',20)
plot([x3,M(1)],[y3,M(2)])
plot([X2,M(1)],[Y2,M(2)])

p1=plot(x4,y4,'r.','markersize',20)
p2=plot(M(1),M(2),'m.','markersize',20)
hold off
grid
axis equal
xlabel('$x$','interpreter','latex','fontsize',14)
ylabel('$y\qquad $ ','interpreter','latex','fontsize',14,'rotation',0)
legend([p1,p2],{'$P_4 = Q_1$',['$P_3+(P_3- P_2)$'  newline  '$ = Q_3 + (Q_2- Q_3)$']},'interpreter','latex','fontsize',10)

axis equal



%% derivatives

dx=(-3+6*t-3*t.^2)*x1+(3-12*t+9*t.^2)*x2+(6*t-9*t.^2)*x3+3*t.^2*x4; %x'(t)
dy=(-3+6*t-3*t.^2)*y1+(3-12*t+9*t.^2)*y2+(6*t-9*t.^2)*y3+3*t.^2*y4; %y'(t)

dX=(-3+6*t-3*t.^2)*X1+(3-12*t+9*t.^2)*X2+(6*t-9*t.^2)*X3+3*t.^2*X4; %X'(t)
dY=(-3+6*t-3*t.^2)*Y1+(3-12*t+9*t.^2)*Y2+(6*t-9*t.^2)*Y3+3*t.^2*Y4; %Y'(t)


d2x=(6-6*t)*x1+(-12+18*t)*x2+(6-18*t)*x3+6*t*x4; %x''(t)
d2y=(6-6*t)*y1+(-12+18*t)*y2+(6-18*t)*y3+6*t*y4; %y''(t)

d2X=(6-6*t)*X1+(-12+18*t)*X2+(6-18*t)*X3+6*t*X4; %X''(t)
d2Y=(6-6*t)*Y1+(-12+18*t)*Y2+(6-18*t)*Y3+6*t*Y4; %Y''(t)

figure(3)
subplot(2,1,1)
plot(t,dx,'b',t+1,dX,'r','linewidth',2)
grid
xlabel('t')
set(gca,'xtick',0:0.5:2)
set(gca,'xticklabel',{'0','0.5','1/0','0.5','1'})
legend({'$x^\prime (t)$','$X^\prime(t)$'},'interpreter','latex','fontsize',11)
title('derivatives')
subplot(2,1,2)
plot(t,dy,'b',t+1,dY,'r','linewidth',2)
xlabel('t')
grid
xlabel('t')
set(gca,'xtick',0:0.5:2)
set(gca,'xticklabel',{'0','0.5','1/0','0.5','1'})
legend({'$y^\prime (t)$','$Y^\prime(t)$'},'interpreter','latex','fontsize',11)

figure(4)
subplot(2,1,1)
plot(t,d2x,'b',t+1,d2X,'r','linewidth',2)
grid
xlabel('t')
set(gca,'xtick',0:0.5:2)
set(gca,'xticklabel',{'0','0.5','1/0','0.5','1'})
legend({'$x^{\prime\prime} (t)$','$X^{\prime\prime}(t)$'},'interpreter','latex','fontsize',11)
title('second derivatives')
subplot(2,1,2)
plot(t,d2y,'b',t+1,d2Y,'r','linewidth',2)
xlabel('t')
grid
xlabel('t')
set(gca,'xtick',0:0.5:2)
set(gca,'xticklabel',{'0','0.5','1/0','0.5','1'})
legend({'$y^{\prime\prime} (t)$','$Y^{\prime\prime}(t)$'},'interpreter','latex','fontsize',11)

