clear
close all

%GIVEN

k = 50
V = 1000
C = 2
D_fi = 2
w = 0.2

la = V / k
% p(t) = C sin(wt) + D

% x'(t) = k * p(t) - k * x(t) / V
% x'(t) + k/V * x(t) = k * C sin(w * t) + k * D
%first degree equation
%a * y'(t) + b * y(t) = A * sin(w * t + theta) + B,y(0)=y0
%solution_formulas.pdf, 1 b ii) p.2

a = 1
b = k / V
A = k * C
B = k * D_fi
y0=0 %q(0)

la=-b/a
K=1/sqrt(b^2+(a*w)^2)
phi=atan2(-a*w,b)
D=y0-K*A*sin(phi)-B/b


% Solution
p = @(t) C * sin(w * t) + D_fi

q=@(t) K*A*sin(w*t+phi)+B/b+D*exp(la*t)

dq=@(t) K*A*cos(w*t+phi)*w+D*exp(la*t)*la %q'(t)=i(t)

tau=-1/la %time constant
tmax=10*tau

t=0:tau/100:tmax;

figure(1)

hold
plot(t, p(t),'blue', 'linewidth', 1.5)
plot(t, q(t) / V,'red','linewidth', 1.5)
grid
xlabel('time (min)')
ylabel('salt concentrations(kg/litra)')
title(['k = ', num2str(k),', V = ', num2str(V), ', C = ', num2str(C), ', D = ', num2str(D_fi), ', \omega  = ', num2str(w)])
legend({'p(t)','x(t)/V'},'fontsize', 12)
hold off

w = 0:0.1:10
K_t = @(w) 1 ./ sqrt(b .^ 2 + (a .* w) .^ 2)
figure(2)
hold
grid
title(['vahvistus K'])

plot(w, K_t(w) / (V / k), 'blue', 'linewidth', 1.5)
hold off

theta_t = @(w) atan2(-a .* w, b)
figure(3)
grid
hold
title(['vahvistus \theta'])
xlabel('kulmataajuus \omega')
plot(w, theta_t(w), 'blue', 'linewidth', 1.5)
hold off

