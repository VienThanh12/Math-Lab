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


