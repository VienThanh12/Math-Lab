%% p.1-6
clear
format long %shows 15 decimals
%format short %4 decimals
N=100000%number of time intervals
t1=1 %start time
t2=2 %end time
dt=(t2-t1)/N %time step

v=@(t) t.^3-10*t.^2 +25*t;
tk=(t1+dt/2):dt:t2; %center points of time intervals
vk=v(tk); %velocities
areas=vk*dt; %areas
area=sum(areas)
%%
figure(1)
bar(tk,vk,1,'g') %1 = width of bar, between 0...1
hold on
t=0:0.01:5;
plot(t,v(t),'b','linewidth',1.5)
%plot(tk,vk,'r.','markersize',12)
plot([t1,t1],[0,v(t1)],'k')
plot([t2,t2],[0,v(t2)],'k')
hold off
xlabel('time t')
ylabel('velocity v(t)')
title(['t_1 = ',num2str(t1),', t_2 = ',num2str(t2),...
    ', N = ',num2str(N),', sum of areas = ',num2str(area)])
xlim([0,5])
grid
%%
%position
s=@(t) 1/4*t.^4-10/3*t.^3+25/2*t.^2;
%distance traveled
s(t2)-s(t1)

%% p.10-12
clear
N=10000
dx=1/N
xk=dx/2:dx:1;
f=xk.^2
sum(f*dx)

%% p.20
clear
f=@(x) 2*x.^2-x-1
F=@(x) 2/3*x^3-1/2*x^2-x
F(2)-F(-1)
%% numerically
N=100
a=-1
b=2
dx=(b-a)/N
x=(a+dx/2):dx:b;%center points
x=a:dx:(b-dx);%left endpoints
x=(a+dx):dx:b;%right endpoints
sum(f(x)*dx)
%% p.21 area of a circle
clear
r=1
a=-r
b=r
f=@(x) sqrt(r^2-x.^2)
N=1000
dx=(b-a)/N
xk=(a+dx/2):dx:b;
sum(f(xk)*dx)
1/2*pi*r^2
%% p.38, volume of a cone
clear
r=2
h=5
N=100 %number of discs
dx=h/N %width of discs
x=dx/2:dx:h; %middle points
rx=r/h*x; %ratio of the cross section
Ax=pi*rx.^2; %area of the cross section
dV=Ax*dx; %volumes of the discs
V=sum(dV)
1/3*pi*r^2*h


%% p.51, length of graph y=x^2
clear
close all
a=0
b=2
dx=(b-a)/1000
x=a:dx:b;
y=x.^2;
s=(asinh(2*x)+2*x.*sqrt(4*x.^2+1))/4;
k=501
subplot(1,2,1)
plot(x,y,x(1:k),y(1:k),'r',x(k),y(k),'r.','markersize',20,'linewidth',1.5)
grid
xlabel('x')
ylabel('y = x^2')

subplot(1,2,2)
plot(x,s,x(k),s(k),'r.','linewidth',1.5,'markersize',20)
grid
xlabel('x')
ylabel('k�yr�n pituus s')
title(['x_0 = ',num2str(x(k)),', s = ',num2str(s(k))])

%% numerically

dy=y(2:k)-y(1:k-1);%vertical distances of consecutive points
ds=sqrt(dx^2+dy.^2); %distances of consecutive points
sk=sum(ds) %length of the graph
s(k)
%% p. 55, square wheel bicycle

clear
close all
xL=asinh(1)
x=-xL:xL/100:xL;
f=@(x) -cosh(x)
df=@(x) -sinh(x)
for k=1:length(x)
%k=50
x0=x(k);
f0=f(x0);
df0=df(x0);
s=sinh(x0)+1;
th=atan2(df0,1);
Ax=x0-s*cos(th);
Ay=f0-s*sin(th);
Bx=x0+(2-s)*cos(th);
By=f0+(2-s)*sin(th);
Cx=Bx+2*cos(th+pi/2);
Cy=By+2*sin(th+pi/2);
Dx=Ax+2*cos(th+pi/2);
Dy=Ay+2*sin(th+pi/2);
plot(x,f(x),'linewidth',2)
hold on
plot([Ax,Bx,Cx,Dx,Ax],[Ay,By,Cy,Dy,Ay],'r','linewidth',2)
plot([-1.5,1.5],[0,0],'k')
plot(x0,f0,'b.','markersize',20)
%keskipiste
Px=(Ax+Cx)/2
Py=(Ay+Cy)/2
plot(Px,Py,'r.','markersize',20)
hold off
grid
axis([-2.5,2.5,-2.5,2.5])
axis square
pause(0.01)
end

%% catenary, p.61
clear
close all
x1=1
y1=4
x2=5
y2=6
AB=sqrt((x2-x1)^2+(y2-y1)^2)
L=8 %length of the rope >AB
%s=length of the curve
%y=a*cosh((x-b)/a)+c
%parameter a=aL from the condition s = L
%by calculating s for a = amin ... amax
amin=0.1
amax=5
a=amin:(amax-amin)/1000:amax;
d=(y1 - y2)./(2*a.*sinh((x1 - x2)./(2*a)));
b=1/2*(x1 + x2 - 2*a.*asinh(d));
s=a.*(sinh((x2-b)./a)-sinh((x1-b)./a));

%find the last a such that s>L
aL_plus=a(s>=L);
aL=aL_plus(end);
%corresponding s
sL_plus=s(s>=L);
sL=sL_plus(end);



figure(1)
plot(a,s,'linewidth',1.5)
hold
plot([0,amax],[AB,AB],'k','linewidth',1.5)
plot([0,aL,aL],[sL,sL,0],'g.-','linewidth',1.5,'markersize',20)
hold off
ylim([0,16])
grid
xlabel('parameter a')
ylabel('length s of the curve')
title(['AB = ',num2str(AB),', L = ',num2str(L),', a = ',num2str(aL),', s = ',num2str(sL)])

a=aL
d=(y1 - y2)/(2*a*sinh((x1 - x2)/(2*a)));
b=1/2*(x1 + x2 - 2*a*asinh(d));
c=y1-a*cosh((x1-b)/a)

%% graph
x=x1:(x2-x1)/100:x2;
y=a*cosh((x-b)/a)+c;

figure(2)
plot(x,y,'linewidth',1.5)
hold
p1=plot(x1,y1,'r.','markersize',20)
p2=plot(x2,y2,'g.','markersize',20)
p3=plot(b,a+c,'b.','markersize',20)
hold off
grid
axis equal
legend([p1,p2,p3],{'A','B','[b,a+c]'},'location','northwest')
title(['L = ',num2str(L),': a = ',num2str(a),', b = ',num2str(b),', c = ',num2str(c)])



%% p.77 curvature <-> direction angle
%paraabola y=x^2
clear
close all
a=0
b=2
dx=(b-a)/1000
x=a:dx:b;
y=x.^2;

s=(asinh(2*x)+2*x.*sqrt(4*x.^2+1))/4;

df=2*x; %f'(x)
d2f=2; %f''(x)
kappa=d2f./(sqrt(1+df.^2)).^3; %curvature k(x)
%direction (polar) angles of the tangent [1,f'(x)]
th=atan2(df,1);

figure(1)
plot(x,y,'linewidth',1.5)
hold
k=501
dk=50
ds=s(k+dk)-s(k)
L=1.2
plot([x(k),x(k)+L*cos(th(k))],[y(k),y(k)+L*sin(th(k))],'r','linewidth',1.5)
plot([x(k+dk),x(k+dk)+L*cos(th(k+dk))],[y(k+dk),y(k+dk)+L*sin(th(k+dk))],'g','linewidth',1.5)
p1=plot(x(1:k),y(1:k),'r','linewidth',1.5)
p2=plot(x(k:k+dk),y(k:k+dk),'k','linewidth',1.5)
p3=plot(x(k),y(k),'r.','markersize',15)
p4=plot(x(k+dk),y(k+dk),'g.','markersize',15)
hold off
grid
axis equal
title({['s = ',num2str(s(k)),', \theta(s) = ',num2str(th(k)),...
       ', \kappa(s) = ',num2str(kappa(k))];...
      ['\Delta s = ',num2str(ds),', \theta(s+\Delta s) = ',num2str(th(k+dk))];...
      ['\Delta\theta = ',num2str(th(k+dk)-th(k)),...
       ', \kappa(s)\Delta s = ',num2str(kappa(k)*ds)]})
legend([p1,p2,p3,p4],{'s','\Delta s','x(s),y(s)','x(s+\Delta s),y(s+\Delta s)'},...
       'fontsize',12,'location','northwest')


%%
figure(2)
subplot(2,1,1)
plot(s,th,'linewidth',1.5)
grid
title('polar angle \theta(s)')
subplot(2,1,2)
plot(s,kappa,'linewidth',1.5)
grid
title('curvature \kappa(s)')
xlabel('length of the curve s')

%% s.79- , Euler spiral
clear
%close all
k0=0 %k(0)
kL=1.5 %k(L)
th0=0*pi/180 %th(0)
thL=220*pi/180 %th(L)
L=2*(thL-th0)/(kL+k0)
%initial point of the curve
x0=0 %x(0)
y0=0 %y(0)


ds=L/100 %step
s=0:ds:L;
ks=k0+(kL-k0)/L*s; %curvature
ths=th0+k0*s+(kL-k0)/(2*L)*s.^2; %direction angle
figure(1)
subplot(2,1,1)
plot(s,ks,'linewidth',1.5)
grid
title('curvature \kappa(s)')

subplot(2,1,2)
plot(s,ths,'linewidth',1.5)
grid
title('direction angle \theta(s)')
xlabel('length s of the cruve ')

%% points on the curve
%move length ds to the direction of the angle ths

N=length(s)
x=zeros(1,N);
y=zeros(1,N);

%initial point
x(1)=x0;
y(1)=y0;

for n=1:N-1;
   x(n+1)=x(n)+ds*cos(ths(n));
   y(n+1)=y(n)+ds*sin(ths(n));
end



figure(2)
plot(x,y,'linewidth',1.5)
hold
plot(x(1),y(1),'g.','markersize',20)
plot(x(end),y(end),'m.','markersize',20)


hold off
grid
axis equal
xlabel('x')
ylabel('y','rotation',0)
title([
     '\kappa(0) = ',num2str(ks(1)),...
     ', \kappa(L) = ',num2str(ks(end)),...
     ', \theta(0) = ',num2str(ths(1)),...
     ', \theta(L) = ',num2str(ths(end),3),...
      ', L = ',num2str(L,3),...
           ])

