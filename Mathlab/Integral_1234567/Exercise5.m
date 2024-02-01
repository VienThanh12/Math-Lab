clear
close all
format long


% GIVEN
v_0 = 10
a = pi / 6
h = 5
g = 9.81

% Solve
% - g / (2 *(v0 * cos(a))^2) * (x^2) + tan(a) * x + h

first = - g / (2 * (v_0 * cos(a)) ^ 2)
second = tan(a)
delta = second ^ 2 - 4 * first * h
x1 = (- second + sqrt(delta)) / (2 * first)
x2 = (- second - sqrt(delta)) / (2 * first)


if(x1 > 0)
  k = x1
else
  k = x2
end

% The length s of the picture

s = -(2^(3/2)*cos(a)^4*v_0^4*asinh((sec(a)^2*sin(2*a)*v_0^2-2*sec(a)^2*g*k)/(2*v_0^2))+(2*cos(a)^2*tan(a)*v_0^2-2*g*k)*sqrt((cos(2*a)+1)*v_0^4-2*sin(2*a)*g*k*v_0^2+2*g^2*k^2)+(-2^(3/2)*cos(a)^4*asinh((sec(a)^2*sin(2*a))/2)-2*cos(a)^2*tan(a)*sqrt(cos(2*a)+1))*v_0^4)/(2^(5/2)*cos(a)^2*g*v_0^2)

% Graph
figure(1)
hold on
x = 0:0.1:k
f = @(x) -g ./ (2 .* (v_0 .* cos(a)) .^ 2) .* (x .^ 2) + tan(a) .* x + h
plot(x, f(x), 'linewidth', 1.5)
xlim([-1, k])
ylim([-2, v_0])
plot([0, k], [0, 0], 'black' ,'linewidth', 1)
plot([0, 0], [0, h + 3], 'black' ,'linewidth', 1)
title(['v_0 = ',num2str(v_0),', a = ',num2str(30), ', h = ', num2str(h), ', s = ', num2str(s)])

xlabel('x')
ylabel('y')
grid
hold off
