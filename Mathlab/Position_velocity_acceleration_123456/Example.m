%% p.3-4
clear
close all
%position
s=@(t) 2*t.^4-5*t.^3+6*t
%velocity
v=@(t) 8*t.^3-15*t.^2+6
%acceleration
a=@(t) 24*t.^2-30*t
%values at time t0
t0=0.5
s(t0)
v(t0)
a(t0)

%graphs when t=0...T
T=2
t=0:T/100:T;

subplot(311)
plot(t,s(t),'linewidth',1.5)
grid
title('position')
subplot(312)
plot(t,v(t),'linewidth',1.5)
grid
title('velocity')
subplot(313)
plot(t,a(t),'linewidth',1.5)
grid
title('acceleration')
xlabel('time t')


%% p.5-6, S-curve
%s(t)=At^5+Bt^4+Ct^3+Dt^2+Et+F
%v(t)=5At^4+4Bt^3+3Ct^2+2Dt+E
%a(t)=20At^3+12Bt^2+6Ct+2D

%s(0)=F=0
F=0
%v(0)=E=0
E=0
%a(0)=2D=0
D=0

%s(1)=A+B+C=1
%v(1)=5A+4B+3C=0
%a(1)=20A+12B+6C=0

%linear system of equations
%coefficient matrix
K=[1,1,1
   5,4,3
   20,12,6]
%right hand side
R=[1
   0
   0]
%solution
X=K\R
A=X(1)
B=X(2)
C=X(3)


s=@(t) A*t.^5+B*t.^4+C*t.^3+D*t.^2+E*t+F
v=@(t) 5*A*t.^4+4*B*t.^3+3*C*t.^2+2*D*t+E
a=@(t) 20*A*t.^3+12*B*t.^2+6*C*t+2*D


t=0:0.01:1;


subplot(311)
plot(t,s(t),'linewidth',1.5)
grid
title('position')
subplot(312)
plot(t,v(t),'linewidth',1.5)
grid
title('velocity')
subplot(313)
plot(t,a(t),'linewidth',1.5)
grid
title('acceleration')
xlabel('time t')



%% p.9-10
clear
close all
%position
s=@(t) sqrt((3+2*t).^2+4^2)-5
%velocity
v=@(t) (6+4*t)./sqrt((3+2*t).^2+4^2)
%acceleration
a=@(t) 64./((3+2*t).^2+4^2).^(3/2)
%values at time t0
t0=2
s(t0)
v(t0)
a(t0)

%graphs when t=0...T
%AP=9->3+2*T=sqrt(9^2-4^2)=sqrt(65)->
T=(sqrt(65)-3)/2
t=0:T/100:T;

subplot(311)
plot(t,s(t),'linewidth',1.5)
grid
xlim([0,T])
title('position')
subplot(312)
plot(t,v(t),'linewidth',1.5)
grid
xlim([0,T])
title('velocity')
subplot(313)
plot(t,a(t),'linewidth',1.5)
grid
xlim([0,T])
title('acceleration')
xlabel('time t')


%% p. 11-12
L=9 %length of rope
h=4 %height
xA=7 %position of A
vA=2 %velocity of A
aA=0 %acceleration of A
sB=sqrt(xA^2+h^2)+h-L %position of B
vB=xA*vA/sqrt(xA^2+h^2) %velocity of B
aB=((xA^3+h^2*xA)*aA+h^2*vA^2)/(xA^2+h^2)^(3/2) %acceleration of B
%compare
s(t0)
v(t0)
a(t0)
%% p.13-14
clear
close all
%angle
th=@(t) atan(3*t/5)
%angular velocity
w=@(t) 15./(9*t.^2+25)
%angular acceleration
alfa=@(t) -270*t./(9*t.^2+25).^2
%values at time t0
t0=2
th0=th(t0)
w0=w(t0)
alfa0=alfa(t0)

%graphs when t=0...T
T=15
t=0:T/100:T;

subplot(311)
plot(t,th(t),'linewidth',1.5)
grid
xlim([0,T])
title('angle \theta(t)')
subplot(312)
plot(t,w(t),'linewidth',1.5)
grid
xlim([0,T])
title('angular velocity \omega(t)')
subplot(313)
plot(t,alfa(t),'linewidth',1.5)
grid
xlim([0,T])
title('angular acceleration \alpha(t)')
xlabel('time t')
%% p.15-16

h=5 %height OP
x=6 %position of P
v=3 %velocity of P
a=0 %acceleration of P
th2=atan(x/h) %angle
w2=h*v/(x^2+h^2) %angular velocity
alfa2=h*((x^2+h^2)*a-2*x*v^2)/(x^2+h^2)^2 %angular acceleration
%compare
th0
w0
alfa0
%% p.17 (slider-crank)
clear
close all
r=0.05;
L=0.12;
w=50*2*pi;

T=2*pi/w; %period
dt=T/100; %time step
t=0:dt:T;
n=length(t);

th=w*t;
co=cos(th);
si=sin(th);
nj=sqrt(L^2-(r*si).^2);

s=r*co+nj; %position
v=-w*r*si-w*r^2*si.*co./nj; %velocity
a=-w^2*r*co-w^2*r^2*(co.^2-si.^2)./nj-...
    w^2*r^4*co.^2.*si.^2./nj.^3; %acceleration


figure(1)
subplot(3,1,1)
plot(t,s,'linewidth',1.5)
grid
title('position s(t)')

subplot(3,1,2)
plot(t,v,'linewidth',1.5)
grid
title('velocity v(t)')

subplot(3,1,3)
plot(t,a,'linewidth',1.5)
grid
title('acceleration a(t)')
xlabel('time t')

%% p.19-20, velocity and acceleration numerically without derivatives
dt=T/10000; %time step
tk=0:dt:T;
n=length(tk);

th=w*tk;
co=cos(th);
si=sin(th);
nj=sqrt(L^2-(r*si).^2);

sk=r*co+nj; %position


%average velocities in time intervals
vk=(sk(2:n)-sk(1:n-1))/dt; %n-1 numbers
%  s(2:n)=[s(2),s(3),...,s(n)]
%s(1:n-1)=[s(1),s(2),...,s(n-1)]

%average accelerations
ak=(vk(2:n-1)-vk(1:n-2))/dt; %n-2 numbers

figure(2)
subplot(3,1,1)
plot(t,s,tk,sk,'r.','linewidth',1.5,'markersize',15)
grid
title('position s(t)')
xlim([0,T])

subplot(3,1,2)
plot(t,v,tk(1:n-1),vk,'r.','linewidth',1.5,'markersize',15)
grid
title('velocity v(t)')
xlim([0,T])

subplot(3,1,3)
plot(t,a,tk(1:n-2),ak,'r.','linewidth',1.5,'markersize',15)
grid
title('acceleration a(t)')
xlabel('time t')
xlim([0,T])


%%  p.21-22
clear
close all

OA=0.7
OB=0.2
AC=1.0
CD=0.5
h=0.4
w=20*pi


T=2*pi/w %period
dt=T/1000
t=0:dt:T;
%coordinates
Ax=0
Ay=-OA
alfa=w*t;
%polar angle of OB is alfa-pi/2
Bx=OB*cos(alfa-pi/2);
By=OB*sin(alfa-pi/2);
th=atan2(By-Ay,Bx-Ax); %polar angle of AB and AC
Cx=Ax+AC*cos(th);
Cy=Ay+AC*sin(th);
Dy=h;
CDx=sqrt(CD^2-(Dy-Cy).^2);
Dx=Cx-CDx;
s=Dx;
n=length(t);
%% velocity and acceleration numerically


vk=(s(2:n)-s(1:n-1))/dt; %n-1 numbers
%  s(2:n)=[s(2),s(3),...,s(n)]
%s(1:n-1)=[s(1),s(2),...,s(n-1)]
ak=(vk(2:n-1)-vk(1:n-2))/dt; %n-2 numbers

figure(2)
subplot(3,1,1)
plot(t,s,'linewidth',1.5)
grid
title('position s(t)')

subplot(3,1,2)
plot(t(1:n-1),vk,'linewidth',1.5)
grid
title('velocity v(t)')

subplot(3,1,3)
plot(t(1:n-2),ak,'linewidth',1.5)
grid
title('acceleration a(t)')
xlabel('time t')
%% animation
figure(3)
for k=1:n
plot([Ax,Cx(k),Dx(k)],[Ay,Cy(k),Dy],'b.-','linewidth',1.5,'markersize',15)
hold
plot([0,Bx(k)],[0,By(k)],'k.-','linewidth',1.5,'markersize',15)
hold off
grid
axis([-4,4,-4,4])
% axis square
axis equal
pause(0.01)

end

%% angular velocity and acceleration of AC
wkAC=(th(2:end)-th(1:end-1))/dt;
akAC=(wkAC(2:end)-wkAC(1:end-1))/dt;

figure(4)
subplot(3,1,1)
plot(t,th,'linewidth',1.5)
grid
title('angle \theta(t) (rad)')

subplot(3,1,2)
plot(t(1:n-1),wkAC,'linewidth',1.5)
grid
title('angular velocity \omega(t) = \theta^\prime(t) (rad/sec)')

subplot(3,1,3)
plot(t(1:n-2),akAC,'linewidth',1.5)
grid
title('angular acceleration \alpha(t) = \omega^\prime(t) (rad/sec^2)')
xlabel('time t')
%% p.27-32
clear
close all

x=@(t) 1.5*sin(0.5*t)
y=@(t) cos(1.25*t)
vx=@(t) 0.75*cos(0.5*t)
vy=@(t) -1.25*sin(1.25*t)
v=@(t) sqrt(vx(t).^2+vy(t).^2)
ax=@(t) -0.375*sin(0.5*t)
ay=@(t) -1.5625*cos(1.25*t)
a= @(t) sqrt(ax(t).^2+ay(t).^2)
aT= @(t) (vx(t).*ax(t)+vy(t).*ay(t))./v(t)
aN= @(t) (vx(t).*ay(t)-vy(t).*ax(t))./v(t)



T=4*pi
t=0:T/1000:T;
t0=7.4

figure(1)
plot(x(t),y(t),'linewidth',1.5)
hold
p2=plot([x(t0),x(t0)+ax(t0)],[y(t0),y(t0)+ay(t0)],'g','linewidth',2)
p3=plot([x(t0),x(t0)+aT(t0)*vx(t0)/v(t0)],[y(t0),y(t0)+aT(t0)*vy(t0)/v(t0)],'c','linewidth',3)
p1=plot([x(t0),x(t0)+vx(t0)],[y(t0),y(t0)+vy(t0)],'r','linewidth',2)
p4=plot([x(t0),x(t0)+aN(t0)*(-vy(t0))/v(t0)],[y(t0),y(t0)+aN(t0)*vx(t0)/v(t0)],'m','linewidth',2)
plot(x(t0),y(t0),'b.','markersize',15)
hold off
grid
axis([-2,2,-2,2])
axis square
set(gca,'xtick',-2:0.5:2)
legend([p1,p2,p3,p4],{'v','a','a_T','a_N'},'fontsize',12)
title(['t = ',num2str(t0),', v = ',num2str(v(t0),3),', a = ',num2str(a(t0),3),...
       ', a_T = ',num2str(aT(t0),3),', a_N = ',num2str(aN(t0),3)])


%% animation
figure(2)
for t0=0:T/1000:T
plot(x(t),y(t),'linewidth',1.5)
hold
plot([x(t0),x(t0)+ax(t0)],[y(t0),y(t0)+ay(t0)],'g','linewidth',2)
plot([x(t0),x(t0)+aT(t0)*vx(t0)/v(t0)],[y(t0),y(t0)+aT(t0)*vy(t0)/v(t0)],'c','linewidth',3)
plot([x(t0),x(t0)+vx(t0)],[y(t0),y(t0)+vy(t0)],'r','linewidth',2)
plot([x(t0),x(t0)+aN(t0)*(-vy(t0))/v(t0)],[y(t0),y(t0)+aN(t0)*vx(t0)/v(t0)],'m','linewidth',2)
plot(x(t0),y(t0),'b.','markersize',15)
hold off
grid
axis([-2,2,-2,2])
axis square
pause(0.01)
end
%%  graphs
figure(3)
subplot(2,1,1)
plot(t,v(t),'linewidth',1.5)
grid
title('|v|')
xlim([0,T])
subplot(2,1,2)
plot(t,a(t),'linewidth',1.5)
grid
title('|a|')
xlabel('time t')
xlim([0,T])

figure(4)
subplot(2,1,1)
plot(t,aT(t),'linewidth',1.5)
grid
title('tangential acceleration a_T')
xlim([0,T])
subplot(2,1,2)
plot(t,aN(t),'linewidth',1.5)
grid
title('normal acceleration a_N')
xlabel('time t')
xlim([0,T])

%% p.33-
clear
close all
OM=5;
MP=3;
w1=2*pi;
w2=4*pi;

T=2*pi/w1;
dt=T/1000;
t=0:dt:T;
n=length(t);

alfa=w1*t;
beta=w2*t;
Mx=OM*cos(alfa);
My=OM*sin(alfa);

x=Mx+MP*cos(alfa+beta);
y=My+MP*sin(alfa+beta);

vx=-OM*w1*sin(alfa)-MP*(w1+w2)*sin(alfa+beta);
vy=OM*w1*cos(alfa)+MP*(w1+w2)*cos(alfa+beta);

ax=-OM*w1^2*cos(alfa)-MP*(w1+w2)^2*cos(alfa+beta);
ay=-OM*w1^2*sin(alfa)-MP*(w1+w2)^2*sin(alfa+beta);

v=sqrt(vx.^2+vy.^2);
a=sqrt(ax.^2+ay.^2);
aT=(vx.*ax+vy.*ay)./v;
aTx=aT.*vx./v;
aTy=aT.*vy./v;
aN=(vx.*ay-vy.*ax)./v;
aNx=aN.*(-vy)./v;
aNy=aN.*vx./v;

%%
%picture at time t(k)

k=381
p=2 %max length of vectors v,a,aT,aN
R=OM+MP+1; %limit for axes
amax=max(a);
vmax=max(v);

figure(1)
plot([x(k),x(k)+p*vx(k)/vmax],[y(k),y(k)+p*vy(k)/vmax],'g','linewidth',3)
hold on
plot([x(k),x(k)+p*ax(k)/amax],[y(k),y(k)+p*ay(k)/amax],'k','linewidth',3)
plot([x(k),x(k)+p*aTx(k)/amax],[y(k),y(k)+p*aTy(k)/amax],'c','linewidth',3)
plot([x(k),x(k)+p*aNx(k)/amax],[y(k),y(k)+p*aNy(k)/amax],'m','linewidth',3)
plot(x,y,'linewidth',2)
plot([0,Mx(k)],[0,My(k)],'b.-','linewidth',2,'markersize',15)
plot([Mx(k),x(k)],[My(k),y(k)],'r.-','linewidth',2,'markersize',15)
hold off
axis([-R,R,-R,R])
axis square
grid
title(['|v| = ',num2str(v(k)),', |a| = ',num2str(a(k)),', a_T = ',num2str(aT(k)),', a_N = ',num2str(aN(k))])
legend('v','a','a_T','a_N')

%% graphs of |v| and |a|
figure(2)
subplot(2,1,1)
plot(t,v,'linewidth',2)
grid
title('|v|')
subplot(2,1,2)
plot(t,a,'linewidth',2)
grid
title('|a|')
xlabel('time t')

figure(3)
subplot(2,1,1)
plot(t,aT,'linewidth',2)
grid
title('tangential acceleration a_T')
subplot(2,1,2)
plot(t,aN,'linewidth',2)
grid
title('normal acceleration a_N')
xlabel('time t')


%% animation

figure(1)
for k=1:n
plot([x(k),x(k)+p*vx(k)/vmax],[y(k),y(k)+p*vy(k)/vmax],'g','linewidth',3)
hold on
plot([x(k),x(k)+p*ax(k)/amax],[y(k),y(k)+p*ay(k)/amax],'k','linewidth',3)
plot([x(k),x(k)+p*aTx(k)/amax],[y(k),y(k)+p*aTy(k)/amax],'c','linewidth',3)
plot([x(k),x(k)+p*aNx(k)/amax],[y(k),y(k)+p*aNy(k)/amax],'m','linewidth',3)
plot(x,y,'linewidth',2)
plot([0,Mx(k),x(k)],[0,My(k),y(k)],'r.-','linewidth',2,'markersize',15)
hold off
grid
axis([-R,R,-R,R])
axis square
pause(0.01)
end



%% numerically, without derivates

%average velocities
vx=(x(2:n)-x(1:n-1))/dt; %n-1 numbers
vy=(y(2:n)-y(1:n-1))/dt; %n-1 numbers

%average accelerations
ax=(vx(2:n-1)-vx(1:n-2))/dt; %n-2 numbers
ay=(vy(2:n-1)-vy(1:n-2))/dt; %n-2 numbers

v=sqrt(vx.^2+vy.^2); %n-1 numbers
a=sqrt(ax.^2+ay.^2); %n-2 numbers
aT=(vx(1:n-2).*ax+vy(1:n-2).*ay)./v(1:n-2); %n-2
aTx=aT.*vx(1:n-2)./v(1:n-2); %n-2
aTy=aT.*vy(1:n-2)./v(1:n-2); %n-2
aN=(vx(1:n-2).*ay-vy(1:n-2).*ax)./v(1:n-2); %n-2
aNx=aN.*(-vy(1:n-2))./v(1:n-2); %n-2
aNy=aN.*vx(1:n-2)./v(1:n-2); %n-2



figure(4)
subplot(2,1,1)
plot(t(1:n-1),v,'linewidth',2)
grid
title('speed v')
subplot(2,1,2)
plot(t(1:n-2),a,'linewidth',2)
grid
title('acceleration a')
xlabel('aika t')



figure(5)
subplot(2,1,1)
plot(t(1:n-2),aT,'linewidth',2)
grid
title('tangential acceleration a_T')
ylim([-400,400])
subplot(2,1,2)
plot(t(1:n-2),aN,'linewidth',2)
grid
title('normal acceleration a_N')
xlabel('time t')
%% p.37, four-bar
clear
close all
OA=1.5
AB=3.5
BC=4
OC=5
w=2*pi
T=2*pi/w
dt=T/360;
t=0:dt:T;
alfa=w*t;
Ox=0
Oy=0
Cx=OC
Cy=0
Ax=OA*cos(alfa);
Ay=OA*sin(alfa);
AC=sqrt((Ax-Cx).^2+(Ay-Cy).^2);
beta1=asin(Ay./AC);
beta2=acos((AC.^2+BC^2-AB^2)./(2*AC*BC));
beta=beta1+beta2;
Bx=Cx-BC*cos(beta);
By=Cy+BC*sin(beta);
Px=(Ax+Bx)/2;
Py=(Ay+By)/2;

x=Px;
y=Py;
n=length(t)
%everage velocities
vx=(x(2:n)-x(1:n-1))/dt; %n-1 lukua
vy=(y(2:n)-y(1:n-1))/dt; %n-1 lukua
%average accelerations
ax=(vx(2:n-1)-vx(1:n-2))/dt; %n-2 lukua
ay=(vy(2:n-1)-vy(1:n-2))/dt; %n-2 lukua
v=sqrt(vx.^2+vy.^2); %n-1 lukua
a=sqrt(ax.^2+ay.^2); %n-2 lukua
aT=(vx(1:n-2).*ax+vy(1:n-2).*ay)./v(1:n-2); %n-2 lukua
aTx=aT.*vx(1:n-2)./v(1:n-2); %n-2 lukua
aTy=aT.*vy(1:n-2)./v(1:n-2); %n-2 lukua
aN=(vx(1:n-2).*ay-vy(1:n-2).*ax)./v(1:n-2); %n-2 lukua
aNx=aN.*(-vy(1:n-2))./v(1:n-2); %n-2 lukua
aNy=aN.*vx(1:n-2)./v(1:n-2); %n-2 lukua

vmax=max(v);
amax=max(a);
p=1.0; %max-length of vectors v,a,aT,aN in the picture

figure(1)
k=241
plot([Ox,Ax(k),Bx(k),Cx],[Oy,Ay(k),By(k),Cy],'k.-','linewidth',1.5,'markersize',20)
hold
plot(x,y,'b','linewidth',1.5)
p1=plot([x(k),x(k)+p*vx(k)/vmax],[y(k),y(k)+p*vy(k)/vmax],'r','linewidth',2)
p2=plot([x(k),x(k)+p*ax(k)/amax],[y(k),y(k)+p*ay(k)/amax],'g','linewidth',2)
p3=plot([x(k),x(k)+p*aTx(k)/amax],[y(k),y(k)+p*aTy(k)/amax],'c','linewidth',2)
p4=plot([x(k),x(k)+p*aNx(k)/amax],[y(k),y(k)+p*aNy(k)/amax],'m','linewidth',2)
plot(x(k),y(k),'b.','markersize',20)
hold off
grid
xlim([-1,6])
axis equal
title(['|v| = ',num2str(v(k)),', |a| = ',num2str(a(k)),...
    ', a_T = ',num2str(aT(k)),', a_N = ',num2str(aN(k))])
legend([p1,p2,p3,p4],'v','a','a_T','a_N')

%% animation

p=1.0;

figure(2)
for k=1:n-2
plot([Ox,Ax(k),Bx(k),Cx],[Oy,Ay(k),By(k),Cy],'k.-','linewidth',1.5,'markersize',20)
hold
plot(x,y,'b','linewidth',1.5)
p1=plot([x(k),x(k)+p*vx(k)/vmax],[y(k),y(k)+p*vy(k)/vmax],'r','linewidth',2)
p2=plot([x(k),x(k)+p*ax(k)/amax],[y(k),y(k)+p*ay(k)/amax],'g','linewidth',2)
p3=plot([x(k),x(k)+p*aTx(k)/amax],[y(k),y(k)+p*aTy(k)/amax],'c','linewidth',2)
p4=plot([x(k),x(k)+p*aNx(k)/amax],[y(k),y(k)+p*aNy(k)/amax],'m','linewidth',2)
plot(x(k),y(k),'b.','markersize',20)
hold off
grid
axis([-2,6,-3,5])
axis square
%legend([p1,p2,p3,p4],'v','a','a_T','a_N')
pause(0.01)
end
%%
figure(3)
subplot(2,1,1)
plot(t(1:n-1),v,'linewidth',1.5)
grid
title('speed $v$','interpreter','latex')
subplot(2,1,2)
plot(t(1:n-2),a,'linewidth',1.5)
grid
title('acceleration $a$','interpreter','latex')
xlabel('aika t')



figure(4)
subplot(2,1,1)
plot(t(1:n-2),aT,'linewidth',1.5)
grid
title('tangential acceleration $a_T$','interpreter','latex')
%ylim([-400,400])
subplot(2,1,2)
plot(t(1:n-2),aN,'linewidth',1.5)
grid
title('normal acceleration $a_N$','interpreter','latex')
xlabel('time $t$','interpreter','latex')

%angular velocity and acceleration of beta
wbeta=(beta(2:n)-beta(1:n-1))/dt;
abeta=(wbeta(2:n-1)-wbeta(1:n-2))/dt;

figure(5)
subplot(3,1,1)
plot(t,beta,'linewidth',1.5)
grid
title('angle \beta')
subplot(3,1,2)
plot(t(1:n-1),wbeta,'linewidth',1.5)
grid
title('angular velocity')
subplot(3,1,3)
plot(t(1:n-2),abeta,'linewidth',1.5)
grid
title('angular acceleration')
xlabel('aika t')

%% p.41- , 2D-arm
%a)
clear
close all
OM=5
MP=3
Ax=1
Ay=6
Bx=6
By=4
T=2

dt=T/1000
t=0:dt:T;
n=length(t)

s=6*(t/T).^5-15*(t/T).^4+10*(t/T).^3; %s(t)

%position of P
x=Ax+s*(Bx-Ax);
y=Ay+s*(By-Ay);

%circle, center [x0,y0], radius r
%angular velocity w
% x0=3
% y0=5
% r=1
% w=2*pi
% x=x0+r*cos(w*t);
% y=y0+r*sin(w*t);

%%
%velocity and acceleration of P

ds=(30*(t/T).^4-60*(t/T).^3+30*(t/T).^2)*1/T; %s'(t)
d2s=(120*(t/T).^3-180*(t/T).^2+60*(t/T))*(1/T)^2; %s''(t)

vx=ds*(Bx-Ax);
vy=ds*(By-Ay);
ax=d2s*(Bx-Ax);
ay=d2s*(By-Ay);

v=sqrt(vx.^2+vy.^2);
a=sqrt(ax.^2+ay.^2);


figure(1)
subplot(2,1,1)
plot(t,v,'linewidth',2)
grid
title('speed |v|')

subplot(2,1,2)
plot(t,a,'linewidth',2)
grid
title('length of accleration |a|')
xlabel('time t (sek)')


%% inverse kinematics: x,y -> alfa ja beta

OP=sqrt(x.^2+y.^2);
phi=atan(y./x);
kO=acos((OM^2+OP.^2-MP^2)./(2*OM*OP));
kM=acos((OM^2+MP^2-OP.^2)./(2*OM*MP));
alfa=phi-kO;
beta=pi-kM;

%angular velocities and accelerations numerically

valfa=(alfa(2:n)-alfa(1:n-1))/dt;
vbeta=(beta(2:n)-beta(1:n-1))/dt;
aalfa=(valfa(2:n-1)-valfa(1:n-2))/dt;
abeta=(vbeta(2:n-1)-vbeta(1:n-2))/dt;


figure(2)
subplot(3,1,1)
plot(t,alfa,'r','linewidth',1.5)
hold on
plot(t,beta,'b','linewidth',1.5)
hold off
grid
title('angles \alpha and \beta (rad)')
legend({'\alpha','\beta'},'fontsize',11)


subplot(3,1,2)
plot(t(1:n-1),valfa,'r','linewidth',1.5)
hold on
plot(t(1:n-1),vbeta,'b','linewidth',1.5)
hold off
grid
title('angular velocities (rad/sec)')
legend({'\alpha','\beta'},'fontsize',11)


subplot(3,1,3)
plot(t(1:n-2),aalfa,'r','linewidth',1.5)
hold on
plot(t(1:n-2),abeta,'b','linewidth',1.5)
hold off
grid
title('angular accelerations (rad/sec^2)')
legend({'\alpha','\beta'},'fontsize',11)
xlabel('time t (sec)')



%% animation

%direct kinematics: alfa,beta -> Px,Py

Mx=OM*cos(alfa);
My=OM*sin(alfa);

Px=Mx+MP*cos(alfa+beta);
Py=My+MP*sin(alfa+beta);

R=OM+MP+1 %limit for axis

figure(3)
for k=1:2:n
  plot(x(k),y(k),'gs','markersize',12)
  hold on
  plot([0,Mx(k)],[0,My(k)],'b-o','linewidth',2)
  plot([Mx(k),x(k)],[My(k),y(k)],'r-o','linewidth',2)
  plot(x,y,'k','linewidth',2)
  plot(Px(k),Py(k),'r*','linewidth',2)
  hold off
  grid
  axis([-1,R,-1,R])
  axis square
  pause(0.01)
end


%% b)
clear
close all
OM=5
MP=3
Ax=1
Ay=6
Bx=6
By=4
T=2

dt=T/1000
t=0:dt:T;
n=length(t)

s=6*(t/T).^5-15*(t/T).^4+10*(t/T).^3; %s(t)

x=Ax
y=Ay
OP=sqrt(x.^2+y.^2);
phi=atan(y./x);
kO=acos((OM^2+OP.^2-MP^2)./(2*OM*OP));
kM=acos((OM^2+MP^2-OP.^2)./(2*OM*MP));
alfaA=phi-kO;
betaA=pi-kM;

x=Bx
y=By
OP=sqrt(x.^2+y.^2);
phi=atan(y./x);
kO=acos((OM^2+OP.^2-MP^2)./(2*OM*OP));
kM=acos((OM^2+MP^2-OP.^2)./(2*OM*MP));
alfaB=phi-kO;
betaB=pi-kM;

alfa=alfaA+s*(alfaB-alfaA);
beta=betaA+s*(betaB-betaA);

figure(1)
plot(t,alfa,'r','linewidth',1.5)
hold
plot(t,beta,'b','linewidth',1.5)
hold off
grid
xlabel('time t')
legend({'\alpha','\beta'},'fontsize',12)


Mx=OM*cos(alfa);
My=OM*sin(alfa);

Px=Mx+MP*cos(alfa+beta);
Py=My+MP*sin(alfa+beta);

vx=(Px(2:n)-Px(1:n-1))/dt;
vy=(Py(2:n)-Py(1:n-1))/dt;
v=sqrt(vx.^2+vy.^2);
ax=(vx(2:n-1)-vx(1:n-2))/dt;
ay=(vy(2:n-1)-vy(1:n-2))/dt;
a=sqrt(ax.^2+ay.^2);

figure(2)
subplot(2,1,1)
plot(t(1:n-1),v,'linewidth',1.5)
grid
title('P:n vauhti v')
subplot(2,1,2)
plot(t(1:n-2),a,'linewidth',1.5)
grid
title('P:n kiihtyvyys a')
xlabel('aika t')

R=OM+MP+1 %limit for axis
figure(3)
plot(Ax,Ay,'r.',Bx,By,'g.','markersize',20)
hold on
plot(Px,Py,'k','linewidth',1.5)
plot([-1,R],[0,0],'k')
plot([0,0],[-1,R],'k')
hold off
grid
axis([-1,R,-1,R])
axis square
legend('A','B')
%%
figure(4)
for k=1:2:n
  plot(Ax,Ay,'r.',Bx,By,'g.','markersize',20)
  hold on
  plot([0,Mx(k)],[0,My(k)],'b-o','linewidth',2)
  hold on
  plot([Mx(k),Px(k)],[My(k),Py(k)],'r-o','linewidth',2)
  plot(Px,Py,'k','linewidth',2)
  plot(Px(k),Py(k),'rs','linewidth',2)
  hold off
  grid
  axis([-1,R,-1,R])
  axis square
  pause(0.01)
end


