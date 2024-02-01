%% p.5
clear
close all

T=5
dt=T/1000
t=0:dt:T;
a= 4/5*t; %a(t)
v= 2/5*t.^2; %v(t)
s= 2/15*t.^3; %s(t)



figure(1)
subplot(3,1,1)
plot(t,a,'linewidth',1.5)
grid
title('kiihtyvyys a(t)')
subplot(3,1,2)
plot(t,v,'linewidth',1.5)
title('nopeus v(t)')
grid
subplot(3,1,3)
plot(t,s,'linewidth',1.5)
grid
title('paikka s(t)')
xlabel('aika t')

%%  p.6

clear
T=5
dt=T/1000
t=0:dt:T;
a= 3*t.^2-20*t+25; %a(t)
v= t.^3-10*t.^2+25*t; %v(t)
s= 1/4*t.^4-10/3*t.^3+25/2*t.^2; %s(t)

% numerically, p.7-8
N=length(t)
vk=zeros(1,N);
sk=zeros(1,N);
vk(1)=0; %initial value v(0)
sk(1)=0; %initial value s(0)
for n=1:N-1
   vk(n+1)=vk(n)+a(n)*dt;
   sk(n+1)=sk(n)+vk(n)*dt+1/2*a(n)*dt^2;
end
figure(2)
subplot(3,1,1)
plot(t,a,'linewidth',1.5)
grid
title('acceleration a(t)')
subplot(3,1,2)
plot(t,v,t,vk,':','linewidth',1.5)
title('velocity v(t)')
grid
subplot(3,1,3)
plot(t,s,t,sk,':','linewidth',1.5)
grid
title('position s(t)')
xlabel('aika t')

%% p.9-10
clear
t1=2
t2=6
T=8
dt=T/1000
t=0:dt:T;

a1=5*t;
a2=-5*t+20;
a3=5*t-40;

v1=5/2*t.^2;
v2=-5/2*t.^2+20*t-20;
v3=5/2*t.^2-40*t+160;

s1=5/6*t.^3;
s2=-5/6*t.^3+10*t.^2-20*t+40/3;
s3=5/6*t.^3-20*t.^2+160*t-1040/3;


a=a1.*(t<=t1)+a2.*((t>t1)&(t<t2))+a3.*(t>=t2);
v=v1.*(t<=t1)+v2.*((t>t1)&(t<t2))+v3.*(t>=t2);
s=s1.*(t<=t1)+s2.*((t>t1)&(t<t2))+s3.*(t>=t2);


%numerically

N=length(t)
vk=zeros(1,N);
sk=zeros(1,N);
vk(1)=0;
sk(1)=0;

for n=1:N-1
    vk(n+1)=vk(n)+a(n)*dt;
    sk(n+1)=sk(n)+vk(n)*dt+1/2*a(n)*dt^2;
end
figure(4)
subplot(3,1,1)
plot(t,a,'linewidth',1.5)
grid
title('kiihtyvyys a(t)')
subplot(3,1,2)
plot(t,v,t,vk,':','linewidth',1.5)
title('nopeus v(t)')
grid
subplot(3,1,3)
plot(t,s,t,sk,':','linewidth',1.5)
grid
title('paikka s(t)')
xlabel('aika t')




%% p.15-16 falling with gravity and air resistance
clear
%close all
h=30
m=2
g=9.81
b=0.5 %coefficient for air resistance


%1) air resistance F(t)=0
dt=0.01 %time step
s1(1)=0 %s(0)
v1(1)=0 %v(0)
k=1
while s1(k)<=h
      F1=0; %air resistance
      a1(k)=1/m*(m*g-F1);%acceleration
      s1(k+1)=s1(k)+v1(k)*dt+1/2*a1(k)*dt^2;
      v1(k+1)=v1(k)+a1(k)*dt;
      k=k+1;
end

t1=0:dt:(k-1)*dt;
T1=t1(end)


%2)  F(t)=b*v
s2(1)=0 %s(0)
v2(1)=0 %v(0)
k=1
while s2(k)<=h
      F2=b*v2(k);
      a2(k)=1/m*(m*g-F2);
      s2(k+1)=s2(k)+v2(k)*dt+1/2*a2(k)*dt^2;
      v2(k+1)=v2(k)+a2(k)*dt;
      k=k+1;
end
t2=0:dt:(k-1)*dt;
T2=t2(end)

%3) F(t)=b*v^2
s3(1)=0 %s(0)
v3(1)=0 %v(0)
k=1
while s3(k)<=h
      F3=b*v3(k)^2;
      a3(k)=1/m*(m*g-F3);
      s3(k+1)=s3(k)+v3(k)*dt+1/2*a3(k)*dt^2;
      v3(k+1)=v3(k)+a3(k)*dt;
      k=k+1;
end
t3=0:dt:(k-1)*dt;
T3=t3(end)

figure(1)
subplot(3,1,1)
plot(t1(1:end-1),a1,'r','linewidth',1.5)
hold on
plot(t2(1:end-1),a2,'g','linewidth',1.5)
plot(t3(1:end-1),a3,'b','linewidth',1.5)
hold off
legend({'F = 0','F = bv','F = bv^2'},'fontsize',10)
title(['m = ',num2str(m),', b = ',num2str(b),', h = ',num2str(h)])
grid
ylabel('acceleration a(t)')

subplot(3,1,2)
plot(t1,v1,'r','linewidth',1.5)
hold on
plot(t2,v2,'g','linewidth',1.5)
plot(t3,v3,'b','linewidth',1.5)
hold off
grid
ylabel('speed v(t)')

subplot(3,1,3)
plot(t1,s1,'r','linewidth',1.5)
hold on
plot(t2,s2,'g','linewidth',1.5)
plot(t3,s3,'b','linewidth',1.5)
hold off
grid
ylabel('distance s(t)')
xlabel('time t')

%% animation
figure(2)
N1=length(t1)
N2=length(t2)
N3=length(t3)

Nmax=max([N1,N2,N3])
s1=[s1,h*ones(1,Nmax-N1)];
s2=[s2,h*ones(1,Nmax-N2)];
s3=[s3,h*ones(1,Nmax-N3)];
v1=[v1,v1(end)*ones(1,Nmax-N1)];
v2=[v2,v2(end)*ones(1,Nmax-N2)];
v3=[v3,v3(end)*ones(1,Nmax-N3)];
t=0:dt:Nmax*dt;

for n=1:Nmax
plot([0,0,4,4],[h,0,0,-1],'k')
hold on
plot(1,h-s1(n),'r.','markersize',20)
plot(2,h-s2(n),'g.','markersize',20)
plot(3,h-s3(n),'b.','markersize',20)
hold off
%ylim([-1,h])
grid

title(['v1 = ',num2str(v1(n),3),', v2 = ',num2str(v2(n),3),...
           ', v3 = ',num2str(v3(n),3),', t = ',num2str(t(n),3)])

pause(0.03)
end
%% p.17- , falling along a curve
clear
close all

%Point P
px=5
py=-4

g=9.81
%% curve 1),  seqment OP
%y=a*x, a=py/px
a=py/px

x1=[0,px];
y1=a*x1;
N1=length(x1)

%fallng time [0,0]->[px,py]
vk=0
t=0

for k=1:N1-1
    %segment P_k->P_(k+1)
    dxk=x1(k+1)-x1(k);
    dyk=y1(k+1)-y1(k);
    sk=sqrt(dxk^2+dyk^2);
    thetak=atan2(dyk,dxk);
    ak=-g*sin(thetak);
    tk=(-vk+sqrt(vk^2+2*ak*sk))/ak;
    t=t+tk;
    vk=vk+ak*tk;
end
T1=t

plot(x1,y1,px,py,'r.','linewidth',1.5,'markersize',20)
grid
axis equal
title(['line y = ax, a = ',num2str(a),', falling time = ',num2str(T1)])

%NOTE:
%aT=-g*sin(thetak)=-g*py/sqrt(px^2+py^2)
s=sqrt(px^2+py^2)
aT=-g*py/s
%s=1/2*aT*T1^2->T1=sqrt(2*s/aT)
sqrt(2*s/aT)

%% curve 2), parabola
%y=ax^2+bx
%py=a*px^2+b*px->a=(py-b*px)/px^2
b=-2 %y'(0)
a=(py-b*px)/px^2

x2=0:px/100:px;
y2=a*x2.^2+b*x2;
N2=length(x2)

%falling time [0,0]->[px,py]
vk=0
t=0

for k=1:N2-1
    %segment P_k->P_(k+1)
    dxk=x2(k+1)-x2(k);
    dyk=y2(k+1)-y2(k);
    sk=sqrt(dxk^2+dyk^2);
    thetak=atan2(dyk,dxk);
    ak=-g*sin(thetak);
    tk=(-vk+sqrt(vk^2+2*ak*sk))/ak;
    t=t+tk;
    vk=vk+ak*tk;
end
T2=t

plot(x2,y2,px,py,'r.','linewidth',1.5,'markersize',20)
grid
axis equal
title(['paraabola y = ax^2+bx,  b = y^\prime(0) = ',num2str(b),', falling time = ',num2str(T2)])
%% curve 3, circle arc, center P0=[x0,0]

alfa=-atan2(py,px) %angle between x-axis and OP
op=sqrt(px^2+py^2)
r=(op/2)/cos(alfa) %radius
x0=r
y0=0
beta=pi-2*alfa; %angle O_P0_P
kk=-pi:beta/100:(-pi+beta); %angles for the arc
%points on the arc
x3=x0+r*cos(kk);
y3=y0+r*sin(kk);
N3=length(x3)

%falling time [0,0]->[px,py]

vk=0
t=0

for k=1:N3-1
    %segment P_k->P_(k+1)
    dxk=x3(k+1)-x3(k);
    dyk=y3(k+1)-y3(k);
    sk=sqrt(dxk^2+dyk^2);
    thetak=atan2(dyk,dxk);
    ak=-g*sin(thetak);
    tk=(-vk+sqrt(vk^2+2*ak*sk))/ak;
    t=t+tk;
    vk=vk+ak*tk;
end
T3=t

plot(x3,y3,px,py,'r.','linewidth',1.5,'markersize',20)
hold
plot([0,x0,px,0],[0,y0,py,0],'k')
plot([x0,px/2],[y0,py/2],'k')
plot(x0,y0,'k.','markersize',20)
kalfa=0:alfa/100:alfa;
p1=plot(r/10*cos(-kalfa),r/10*sin(-kalfa),'linewidth',2)
p2=plot(x0+r/10*cos(kk),y0+r/10*sin(kk),'linewidth',2)
hold off
grid
axis equal
title(['circle arc, falling time = ',num2str(T3)])
legend([p1,p2],{'\alpha','\beta'},'fontsize',13)
%% cycloid through point P=[px,py]

%x = r*(th-sin(th))
%y = r*(-1+cos(th))

%find radius r and angle th0 such that
%px = r*(th0-sin(th0))
%py = r*(-1+cos(th0))
%i.e
%py/px=(-1+cos(th0))/(th0-sin(th0))
%i.e
%py*((th0-sin(th0))-px*(-1+cos(th0))=0
%find th0 by calculating the values of the left hand side
%py*((th-sin(th))-px*(-1+cos(th))
th=0:pi/1000:2*pi;

lhs=py*(th-sin(th))-px*(-1+cos(th));
thplus=th(lhs>=0); %angles for which lhs>0
thminus=th(lhs<0); %angles for which lhs<0
th0=(thplus(end)+thminus(1))/2


figure(1)
plot(th,lhs,'b','linewidth',2)
hold
plot([0,2*pi],[0,0],'k')
plot(th0,0,'r.','markersize',20)
hold off
grid
xlim([0,2*pi])
xlabel('\theta','fontsize',14)
title('left hand side of the equation')

%solve r from px = r0*(th0-sin(th0))
r=px/(th0-sin(th0))
%%
th=0:th0/100:th0;
x4=r*(th-sin(th));
y4=r*(-1+cos(th));

N4=length(x4)

%falling time [0,0]->[px,py]
vk=0
t=0

for k=1:N4-1
    %segment P_k->P_(k+1)
    dxk=x4(k+1)-x4(k);
    dyk=y4(k+1)-y4(k);
    sk=sqrt(dxk^2+dyk^2);
    thetak=atan2(dyk,dxk);
    ak=-g*sin(thetak);
    tk=(-vk+sqrt(vk^2+2*ak*sk))/ak;
    t=t+tk;
    vk=vk+ak*tk;
end
T4=t

figure(2)
plot(x4,y4,'b','linewidth',1.5)
hold
k=51
kpx=r*th(k);
kpy=-r;
plot(kpx+r*cosd(0:360),kpy+r*sind(0:360),'k','linewidth',1.5)
plot(x4(k),y4(k),'b.','markersize',20)
plot([kpx,kpx,x4(k)],[0,kpy,y4(k)],'k')
plot(kpx,kpy,'k.','markersize',20)
plot(px,py,'r.','markersize',20)
plot([0,px+r],[0,0],'k')
hold off
grid
axis([-1,8,-7,2])
axis square
title(['r = ',num2str(r),', \theta_0 = ',num2str(th0),', falling time = ',num2str(T4)])
%% animation, calculate coordinates of the falling point with time step dt
dt=0.00001

%curve 1
vt=0 %speed
k=1 %segment number
m=1 %time number
t=0 %time value
M=10^6;
%coordinates
xt1=zeros(1,M);
yt1=zeros(1,M);
xt1(1)=0;
yt1(1)=0;
while xt1(m)<px
   %segment k
   dxk=x1(k+1)-x1(k);
   dyk=y1(k+1)-y1(k);
   sk=sqrt(dxk^2+dyk^2);
   thetak=atan2(dyk,dxk);
   ak=-g*sin(thetak);
   stk=0; %distance travelled on segment k
   while stk<=sk; %on segment k
   dstk=vt*dt+1/2*ak*dt^2; %distance trabvelld during time step
   dxtk=dstk*cos(thetak); %horizontal
   dytk=dstk*sin(thetak); %vertical
   xt1(m+1)=xt1(m)+dxtk; %new x-coordinate
   yt1(m+1)=yt1(m)+dytk; %new y-coordinate
   stk=stk+dstk;
   m=m+1;
   vt=vt+ak*dt;
   t=t+dt;
   end
   k=k+1;
end
m1=m
t1=t
%drop the extra zeros from the end
xt1=xt1(1:m1);
yt1=yt1(1:m1);

%NOTE:
%aT=-g*sin(theta)=-g*py/sqrt(px^2+py^2)
%op=sqrt(px^2+py^2)
%aT=-g*py/op
%st=1/2*aT*t.^2;
%xt=st.*cos(theta)=st.*px/op=-1/2*g*px*py/op^2*t.^2;
%yt=st.*sin(theta)=st.*py/op=-1/2*g*py^2/op^2*t.^2;


%curve 2
vt=0 %speed
k=1 %segment number
m=1 %time number
t=0 %time value
M=10^6;
%coordinates
xt2=zeros(1,M);
yt2=zeros(1,M);
xt2(1)=0;
yt2(1)=0;
while xt2(m)<px
   %segment k
   dxk=x2(k+1)-x2(k);
   dyk=y2(k+1)-y2(k);
   sk=sqrt(dxk^2+dyk^2);
   thetak=atan2(dyk,dxk);
   ak=-g*sin(thetak);
   stk=0; %distance travelled on segment k
   while stk<=sk; %on segment k
   dstk=vt*dt+1/2*ak*dt^2; %distance trabvelld during time step
   dxtk=dstk*cos(thetak); %horizontal
   dytk=dstk*sin(thetak); %vertical
   xt2(m+1)=xt2(m)+dxtk; %new x-coordinate
   yt2(m+1)=yt2(m)+dytk; %new y-coordinate
   stk=stk+dstk;
   m=m+1;
   vt=vt+ak*dt;
   t=t+dt;
   end
   k=k+1;
end
m2=m
t2=t
%drop the extra zeros from the end
xt2=xt2(1:m2);
yt2=yt2(1:m2);

%curve 3
vt=0 %speed
k=1 %segment number
m=1 %time number
t=0 %time value
M=10^6;
%coordinates
xt3=zeros(1,M);
yt3=zeros(1,M);
xt3(1)=0;
yt3(1)=0;
while xt3(m)<px
   %segment k
   dxk=x3(k+1)-x3(k);
   dyk=y3(k+1)-y3(k);
   sk=sqrt(dxk^2+dyk^2);
   thetak=atan2(dyk,dxk);
   ak=-g*sin(thetak);
   stk=0; %distance travelled on segment k
   while stk<=sk; %on segment k
   dstk=vt*dt+1/2*ak*dt^2; %distance trabvelld during time step
   dxtk=dstk*cos(thetak); %horizontal
   dytk=dstk*sin(thetak); %vertical
   xt3(m+1)=xt3(m)+dxtk; %new x-coordinate
   yt3(m+1)=yt3(m)+dytk; %new y-coordinate
   stk=stk+dstk;
   m=m+1;
   vt=vt+ak*dt;
   t=t+dt;
   end
   k=k+1;
end
m3=m
t3=t
%drop the extra zeros from the end
xt3=xt3(1:m3);
yt3=yt3(1:m3);

%curve 4
vt=0 %speed
k=1 %segment number
m=1 %time number
t=0 %time value
M=10^6;
%coordinates
xt4=zeros(1,M);
yt4=zeros(1,M);
xt4(1)=0;
yt4(1)=0;
while xt4(m)<px
   %segment k
   dxk=x4(k+1)-x4(k);
   dyk=y4(k+1)-y4(k);
   sk=sqrt(dxk^2+dyk^2);
   thetak=atan2(dyk,dxk);
   ak=-g*sin(thetak);
   stk=0; %distance travelled on segment k
   while stk<=sk; %on segment k
   dstk=vt*dt+1/2*ak*dt^2; %distance trabvelld during time step
   dxtk=dstk*cos(thetak); %horizontal
   dytk=dstk*sin(thetak); %vertical
   xt4(m+1)=xt4(m)+dxtk; %new x-coordinate
   yt4(m+1)=yt4(m)+dytk; %new y-coordinate
   stk=stk+dstk;
   m=m+1;
   vt=vt+ak*dt;
   t=t+dt;
   end
   k=k+1;
end
m4=m
t4=t
%drop the extra zeros from the end
xt4=xt4(1:m4);
yt4=yt4(1:m4);
%% tehd��n kaikista putoamisista yht� pitk�t
mmax=max([m1,m2,m3,m4])
xt1=[xt1,xt1(end)*ones(1,mmax-m1)];
yt1=[yt1,yt1(end)*ones(1,mmax-m1)];
xt2=[xt2,xt2(end)*ones(1,mmax-m2)];
yt2=[yt2,yt2(end)*ones(1,mmax-m2)];
xt3=[xt3,xt3(end)*ones(1,mmax-m3)];
yt3=[yt3,yt3(end)*ones(1,mmax-m3)];
xt4=[xt4,xt4(end)*ones(1,mmax-m4)];
yt4=[yt4,yt4(end)*ones(1,mmax-m4)];

%% animation
figure(3)
for m=1:500:mmax
plot(x1,y1,'r',x2,y2,'g',x3,y3,'b',x4,y4,'k','linewidth',1.5,'markersize',20)
hold on
plot(xt1(m),yt1(m),'r.','markersize',20)
plot(xt2(m),yt2(m),'g.','markersize',20)
plot(xt3(m),yt3(m),'b.','markersize',20)
plot(xt4(m),yt4(m),'k.','markersize',20)
hold off
grid
axis equal
legend('1','2','3','4')
pause(0.01)
end
title(['t_1 = ',num2str(t1),', t_2 = ',num2str(t2),', t_3 = ',num2str(t3),', t_4 = ',num2str(t4)])


%% animation in MATLAB, method 2
figure(3)
plot(x1,y1,'r',x2,y2,'g',x3,y3,'b',x4,y4,'k','linewidth',1.5,'markersize',20)
hold on
p1=plot(xt1(1),yt1(1),'r.','markersize',20)
p2=plot(xt2(1),yt2(1),'g.','markersize',20)
p3=plot(xt3(1),yt3(1),'b.','markersize',20)
p4=plot(xt4(1),yt4(1),'k.','markersize',20)
hold off
grid
axis equal
legend('1','2','3','4')

for m=1:2:mmax
    p1.XData=xt1(m);
    p1.YData=yt1(m);
    p2.XData=xt2(m);
    p2.YData=yt2(m);
    p3.XData=xt3(m);
    p3.YData=yt3(m);
    p4.XData=xt4(m);
    p4.YData=yt4(m);
    drawnow limitrate;
end
title(['t_1 = ',num2str(t1),', t_2 = ',num2str(t2),', t_3 = ',num2str(t3),', t_4 = ',num2str(t4)])


%% p.21-  tractrix
clear
V=1 %speed of front wheel
L=1.5 %length of the bike
T=10 %end time
dt=T/1000 %time step
t=0:dt:T;
N=length(t);


%trajectory of rear wheel
x=zeros(1,N);
y=zeros(1,N);

%1) front wheel moves to the right on x-axis with speed V
X=V*t; %X(t)
Y=zeros(1,N); %Y(t)=0
VX=V*ones(1,N); %X'(t)=V
VY=zeros(1,N); %Y'(t)=0
%starting point of the rear wheel
x(1)=0;
y(1)=L;

% %2) front wheel moves on a circle with speed V
% r=2
% X=r*cos(V*t/r);
% Y=r*sin(V*t/r);
% VX=-V*sin(V*t/r);%X'(t)
% VY=V*cos(V*t/r);%Y'(t)
% %starting point of the rear wheel
% x(1)=r;
% y(1)=-L;


s=0 %distance travelled by the rear wheel
v=zeros(1,N-1); %seed of rear wheel
for k=1:N-1
   vx(k)=((X(k)-x(k))*VX(k)+(Y(k)-y(k))*VY(k))/L^2*(X(k)-x(k)); %x'(t)
   vy(k)=((X(k)-x(k))*VX(k)+(Y(k)-y(k))*VY(k))/L^2*(Y(k)-y(k)); %y'(t)
   dx=vx(k)*dt;
   dy=vy(k)*dt;
   x(k+1)=x(k)+dx;
   y(k+1)=y(k)+dy;
   s=s+sqrt(dx^2+dy^2);
   v(k)=sqrt(vx(k)^2+vy(k)^2);
end

%% animation

figure(1)
for k=1:N
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
%% speed of rear wheel
figure(2)
plot(t(1:N-1),v,'linewidth',1.5)
grid
xlabel('aika t')
title('takarenkaan vauhti v')








