clear
close all
format long


% t = 0 . . . 2π

%Given
R = 5
L = 3


% Derivative of x(t)
x = @(t) -3*L*sin(3*t)-R*sin(t)
% Derivative of y(t)
y = @(t) 3*L*cos(3*t)+R*cos(t)



N = 1000%number of time intervals
t1 = 0 %start time
t2 = 2*pi %end time
dt = (t2-t1)/N %time step

f = @(t) sqrt((x(t) .^ 2) + (y(t) .^ 2))

tk = (t1+dt/2):dt:t2; %center points of time intervals
vk = f(tk); %velocities
areas = vk*dt; %areas
area = sum(areas)

% graph

figure(1)
hold on
axis equal
t = 0:0.01:(2*pi)
u = @(t) R * cos(t) + L * cos(3 * t)
v = @(t) R * sin(t) + L * sin(3 * t)
plot(u(t),v(t),'blue' ,'linewidth', 1.5)
xlabel('x')
ylabel('y')
title(['R = ', num2str(R), ', L = ', num2str(L), ', s = ', num2str(area)])

grid
hold off
