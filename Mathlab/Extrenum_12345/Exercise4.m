clear
close all

% PQ = sqrt(h ^ 2 + x ^ 2)
% QB = L - x

% t(x)  = sqrt(h ^ 2 + x ^ 2) + (L - x) / v

L = 10
h = 4
v = 1.05

t = @(x) sqrt(h .^ 2 + x .^ 2) + (L - x) ./ v
% t' = x/sqrt(x^2+h^2)-1/v
% t' = 0
% x = sqrt(h ^ 2 / (v ^ 2 - 1))
w = sqrt(h ^ 2 / (v ^ 2 - 1))
t(w)

if (w > L)
  w = L
end

figure(1)
hold on

x = 0:0.01:L

%[tmin,indmin] = min(t(x)) %min value of vector f is its element number indmin
%xmin = x(indmin)

plot(x, t(x), 'linewidth', 1.5)
plot(w, t(w),'r.','markersize',15)
ylabel('aika t reitillä PQB')
xlabel('AQ')
title(['L = ', num2str(L), ', h = ', num2str(h), ', v = ', num2str(v), ', min t(', num2str(w), ') = ', num2str(t(w))])
grid
hold off

figure(2)
hold on
plot(0, h, 'b.' ,'markersize', 15)
plot(w, 0, 'b.' ,'markersize', 15)

plot([0, w], [h, 0], 'b','linewidth',1.7)
plot([L, w], [0, 0], 'b','linewidth',1.7)
plot([0, w], [0, 0], 'black','linewidth',1)
title(['nopein retti, AQ = ', num2str(w)])
ylim([-1,5])
grid
hold off
