clear
close all

#GIVEN

R = 1
L = 3
U = 5

% First degree equation
# L . i'(t) + R . i(t) = u(t)
%a . y'(t) + b . y(t) = A, y(0) = y0
% U_R = R * i(t) and U_L = L * i'(t)
a = L
b = R
A = U

y0 = 0
la= -b / a

%solution
y = @(t) A / b + (y0 - A / b) * exp(la * t)
y_d = @(t) (y0 - A / b) * exp(la * t) * la

tau = -1 / la %time constant
tmax = 5 * tau
t=0:tau/100:tmax;

U_R = @(t) R * y(t)
U_L = @(t) L * y_d(t)

plot(t, U_R(t),'linewidth',1.5)
hold on
plot(t, U_L(t),'linewidth',1.5)
hold off
grid
xlabel('time t')
xlim([0, tmax])
title(['R = ', num2str(R),', L = ', num2str(L), ', U = ', num2str(U)])
legend({'U_R','U_L'},'fontsize',12)

