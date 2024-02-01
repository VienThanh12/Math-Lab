clear
close all

format long

a = 0.7
h = 4

area = (2 * h * sqrt(h / a)) / 3

%Coordinates of the centroid
xcentroid = 1 ./ area .* -((h - 8) * h) / (4 * a)
%y
ycentroid = 1 ./ area .* (2 * h ^ 2 * sqrt(h / a)) / 5

% graph
hold on
plot([0, 0], [0, h], 'black','linewidth', 1.7)
plot([0, sqrt(h ./ a)], [h, h], 'b','linewidth', 1.7)


plot([-2, h], [0, 0], 'black','linewidth', 0.5)
plot([0, 0], [-1, 5], 'black','linewidth', 0.5)


x = 0:0.1:sqrt(h ./ a) + 0.1
f = @(x) (a .* (x .^ 2))
plot(x, f(x), 'r','linewidth', 1.4)
xlim([-2, 4])
ylim([-1, 5])
xlabel('x')
grid
title(['a = ', num2str(a) ,', h = ',num2str(h), ': A = ', num2str(area), ', xp = ', num2str(xcentroid), ', yp = ', num2str(ycentroid)])
axis equal
plot(xcentroid, ycentroid, 'black.' ,'markersize', 17)
hold off


