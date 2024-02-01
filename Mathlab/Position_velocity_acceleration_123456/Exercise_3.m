clear
close all
format long
%GIVEN

w = 2 .* pi
T = 2 .* pi ./ w
dt = T ./ 100
t = 0:dt:T

E_D = 400
O_A = sqrt(100 .^ 2 + 200 .^ 2)
O_D = 600
A_D = sqrt(100 .^ 2 + 400 .^ 2)
A_B = sqrt(300 .^ 2 + 50 .^ 2)
C_x = 200
C_y = 100
B_C = 50

% Solve
alfa = w .* t
% Coordinates of B
B_x = (C_x + sin(alfa) .* B_C)
B_y = (C_y + cos(alfa) .* B_C)

% OB and B
O_B = sqrt(B_x .^ 2 + B_y .^ 2)
B = acos(B_x ./ O_B)

% AOB and AOD
A_O_B = acos(((A_B .* A_B) - (O_A .* O_A) - (O_B .* O_B)) ./ (- 2 .* O_A .* O_B))
A_O_D = (acos(((A_D .* A_D) - (O_A .* O_A) - (O_D .* O_D)) ./ (- 2 .* O_A .* O_D))) % h

% polar angle
polar_angle = (A_O_B - A_O_D + B)

u = sin(polar_angle)
u_c =  cos(polar_angle)

% Final result
D_x = O_D .* u_c
E_x = (D_x + sqrt(E_D .^ 2 - ((O_D .* u - 400) .^ 2))) / 1000

% graph
figure(1)
n = length(w .* t);
vk = (E_x(2:n) - E_x(1:n-1)) ./ dt %n-1 numbers
ak=(vk(2:n-1)-vk(1:n-2))/dt; %n-2 numbers

subplot(3, 1, 1)
plot(t, E_x, 'linewidth', 1.5)
title('E:n paikka s = Ex (m)')
grid
subplot(3, 1, 2)
plot(t(1:n-1),vk, 'linewidth', 1.5)
grid
title('E:n nopeus v (m / s)')
subplot(3,1,3)
plot(t(1:n-2),ak,'linewidth',1.5)
grid
title('E:n kiihtyvyys a (m / s2)')
xlabel('aika t (s)')
hold off

