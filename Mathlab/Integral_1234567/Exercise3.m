clear
close all
format long

% v1 = 500, p1 = 30
v_1 = 500, p_1 = 30
% v2 = 50, p2 = 475.4679577
v_2 = 50, p_2 = 475.4679577
% v3 = 50, p3 = 600
v_3 = 50, p_3 = 600
% v4 = 100, p4 = 600
v_4 = 100, p_4 = 600
% v5 = 200, p5 = 300
v_5 = 200, p_5 = 300
% v6 = 500, p6 = 57.65397283
v_6 = 500, p_6 = 57.65397283
% k1 = 1.2 k2 = 1.8
k_1 = 1.2, k_2 = 1.8



% v2 - v4  | p_3 - p1 * (v1 / v2) ^ k1
sum1 = (v_4^k_1*((k_1-1)*p_3*v_2^k_1*v_4+(1-k_1)*p_3*v_2^(k_1+1)-p_1*v_1^k_1*v_2)+p_1*v_1^k_1*v_2^k_1*v_4)/((k_1-1)*v_2^k_1*v_4^k_1)
% v4 - v5 | p4 * (v4 / v5) - p1 * (v1 / v2) ^ k1
sum2 = ((k_1-1)*p_4*v_4^(k_1+1)*v_5^k_1*log(v_5)+((1-k_1)*p_4*v_4^(k_1+1)*log(v_4)-p_1*v_1^k_1*v_4)*v_5^k_1+p_1*v_1^k_1*v_4^k_1*v_5)/((k_1-1)*v_4^k_1*v_5^k_1)
% v5 - v1 |  p5 * (v5 / v6) ^ k2 - p1 * (v1 / v2) ^ k1

sum3 = -((k_1-1)*p_5*v_1*v_5^(k_2+k_1)+v_5^k_1*((1-k_1)*p_5*v_1^k_2*v_5+(1-k_2)*p_1*v_1^(k_2+1))+(k_2-1)*p_1*v_1^(k_2+k_1)*v_5)/((k_1-1)*(k_2-1)*v_1^k_2*v_5^k_1)

% Area of the drawing
area = sum1 + sum2 + sum3

% Graph

figure(1)
hold on
x = v_5:1:v_6
y = @(x) p_5 .* (v_5 ./ x) .^ (k_2)
plot(x, y(x), 'b', 'linewidth', 1.5)

x_1 = v_2:1:v_6
y_1 = @(x_1) p_1 .* (v_1 ./ x_1) .^ (k_1)
plot(x_1, y_1(x_1), 'b', 'linewidth', 1.5)

x_2 = v_4:1:v_5
y_2 = @(x_2) p_4 .* (v_4 ./ x_2)
plot(x_2, y_2(x_2), 'b', 'linewidth', 1.5)


plot([v_1, v_6], [p_1, p_6], 'b', 'linewidth', 1.5)
plot([v_3, v_4], [p_3, p_4], 'b', 'linewidth', 1.5)
plot([v_2, v_3], [p_2, p_3], 'b', 'linewidth', 1.5)

plot(v_1, p_1, 'r.', 'markersize', 15)
plot(v_2, p_2, 'r.', 'markersize', 15)
plot(v_3, p_3, 'r.', 'markersize', 15)
plot(v_4, p_4, 'r.', 'markersize', 15)
plot(v_5, p_5, 'r.', 'markersize', 15)
plot(v_6, p_6, 'r.', 'markersize', 15)

xlabel('tilavuus V')
ylabel('paine p')
title(['v1 = ',num2str(v_1),', p1 = ',num2str(p_1), ', v2 = ', num2str(v_2), ', p3 = ', num2str(p_3), ', v4 = ', num2str(v_4), ', v5 = ', num2str(v_5), ', k1 = ', num2str(k_1), ', k2 = ', num2str(k_2), ', Ala = ', num2str(area)])
grid
hold off
