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

