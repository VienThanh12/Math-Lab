clear
close all
# GIVEN
x = 3
y = 4
v = 2


z = sqrt(x ^ 2 + y ^ 2)
%position
s = @(t) sqrt((z + 2 * t) .^ 2 - y .^ 2)
%velocity
v = @(t) (2 .* (2 .* t + z)) ./ sqrt((2 .* t + z) .^2 - y .^2)
%acceleration
a = @(t) -(4 .* y .^ 2) ./ ((2 .* t + z) .^ 2 - y .^ 2) .^ (3 ./ 2)
%values at time t0
t0=2
s(t0)
v(t0)
a(t0)


%graphs when t=0...T
T=2
t=0:T/100:T;


subplot(311)
plot(t,s(t),'linewidth',1)
grid
xlim([0,T])
title('position s(t)')
subplot(312)
plot(t,v(t),'linewidth',1)
grid
xlim([0,T])
title('velocity v(t)')
subplot(313)
plot(t,a(t),'linewidth',1)
grid
xlim([0,T])
title('acceleration a(t)')
xlabel('time t')

