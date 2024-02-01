clear
close all

E = 1
t = 15
k = 3
f = @(R) ((E .* R) / (R + k) - (E .* R) / (R + t));

% df = -(E*(t-k)*(R^2-k*t))/((R+k)^2*(R+t)^2)
% R = - sqrt(k * t) or sqrt(k * t)
% R > 0
% R = sqrt(k * t sub to the equation
R = sqrt(k *t) % R max
f(sqrt(k * t)) % f max

x=0:0.01:60;

u = @(x) ((E .* x) ./ (x + k) - ((E .* x) ./ (x + t)))

%[umax,indmax]=max(u(x)) %largest value of vector f is its element number indmin
%xmax=x(indmax)

hold on
plot(x, u(x), 'linewidth',1.5)
grid
plot(R,f(R),'r.','markersize',15);
title(['E = ', num2str(E), ', Rt = ', num2str(t), ', Rk = ', num2str(k), ' , max u(', num2str(R), ') = ', num2str(f(R))])
ylabel('jännite u')
xlabel('R')
hold off
