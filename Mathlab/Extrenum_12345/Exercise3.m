clear
close all

% MA + MB = 2 * sqrt(h ^ 2 + (b/2) ^ 2)
% MC = sqrt(a ^ 2 - (b/2) ^ 2) - h
% MA + MB + MC = 2 * sqrt(h ^ 2 + (b/2) ^ 2) + sqrt(a ^ 2 - (b/2) ^ 2) - h


a = 5
b = 6

f = @(h) 2 .* sqrt(h .^ 2 + (b ./ 2) .^ 2) + sqrt(a .^ 2 - (b ./ 2) .^ 2) - h
% df = (2*h)/sqrt(h^2+b^2/4)-1
% df = 0 => h = sqrt(b ^ 2 / 12)
u = sqrt(b ^ 2 / 12)
f(u)

k = sqrt(a ^ 2 - (b / 2) ^ 2)

figure(1)

hold on
h = 0:0.01:k


%[fmin,indmin] = min(f(h)) %min value of vector f is its element number indmin
%hmin = h(indmin)

if u >= k
    u = k
end

plot(h, f(h), 'linewidth', 1.5)
plot(u, f(u),'r.','markersize',15)
xlabel('h')
ylabel('s')
title(['a = ', num2str(a), ', b = ', num2str(b), ', h0 = ', num2str(u)])
grid
hold off

figure(2)
hold on

plot([-b / 2, b / 2], [0, 0], 'b','linewidth',1.5)

if( u != k)
    plot([-b / 2, 0], [0, sqrt(a ^ 2 - (b / 2) ^ 2)], 'b','linewidth',1.5)
    plot([b / 2, 0], [0, sqrt(a ^ 2 - (b / 2) ^ 2)], 'b','linewidth',1.5)
    plot(0, sqrt(a ^ 2 - (b / 2) ^ 2), 'b.' ,'markersize', 15)
end

plot(0, u, 'r.' ,'markersize', 17)

plot([0, 0], [u, sqrt(a ^ 2 - (b / 2) ^ 2)], 'r','linewidth',1.7)
plot([0, -b/2], [u, 0], 'r','linewidth',1.7)
plot([0, b/2], [u, 0], 'r','linewidth',1.7)
plot(-b/2, 0, 'b.' ,'markersize', 15)
plot(b/2, 0, 'b.' ,'markersize', 15)

grid
hold off
