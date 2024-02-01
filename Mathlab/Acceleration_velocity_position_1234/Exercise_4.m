clear
close all

% GIVEN


V = 1 % speed of target
v = 1.1 % speed of the pursuer

T = 10 %end time
dt = T / 1000 %time step
t = 0:dt:T;
N = length(t);

% SOLVE

% 1) target moves to the right on x-axis with speed V
X = V * t; %X(t) = Vt
Y = zeros(1, N); %Y(t) = 0

VX = V * ones(1, N); %X'(t)=V
VY = zeros(1, N); %Y'(t)=0

%trajectory of the pursuer
x = zeros(1,N) % zero from 1 to N
y = zeros(1,N)

x(1) = 0
y(1) = 1


for k=1:N-1
   vx(k) = v * (X(k) - x(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
   vy(k) = v * (Y(k) - y(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
   dx = vx(k) * dt;
   dy = vy(k) * dt;
   x(k + 1) = x(k) + dx;
   y(k + 1) = y(k) + dy;
end

% 2) target moves on a circle

% GIVEN

%V = 1
%v = 1.1

%r = 1

%trajectory of the target
%X=r*cos(V*t/r);
%Y=r*sin(V*t/r);
%X(1) = r
%Y(1) = 0
%VX=-V*sin(V*t/r);%X'(t)
%VY=V*cos(V*t/r);%Y'(t)

%trajectory of the pursuer

%x = zeros(1,N) % zero from 1 to N
%y = zeros(1,N)

%x(1) = 0
%y(1) = 0

%for k=1:N-1
%   vx(k) = v * (X(k) - x(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
%   vy(k) = v * (Y(k) - y(k)) / (sqrt((X(k) - x(k)) ^ 2 + (Y(k) - y(k)) ^ 2));
%   dx = vx(k) * dt;
%   dy = vy(k) * dt;
%   x(k + 1) = x(k) + dx;
%   y(k + 1) = y(k) + dy;
%end


figure(1)
for k = 1:N
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
