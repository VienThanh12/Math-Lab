clear
close all

#GIVEN

J = 0.1
L = 0.5
R = 0.3
b_fi = 0.08
K = 0.1
Kp = 5
Ki = 2
Kd = 0.001


a = J * L
b = R * J + b_fi * L
c = b_fi * R + K * K
A = K

% a ω′′(t) + b ω′(t) + c ω(t) = K u(t)
% w'' = (K u(t) - c w(t) - b w'(t))  / a
dt = 0.01 %time step
N = 1000
t = 0:dt:N * dt;
r = ones(1, N + 1);%setpoint
y = zeros(1, N + 1);
u = zeros(1, N);
e = zeros(1, N);
y(1) = 0
dy(1) = 0
integral = 0
previous_error = 0


for k = 1:N
    error=r(k)-y(k);
    e(k)=error;
    integral=integral+error*dt;
    derivative=(error-previous_error)/dt;
    u(k)=Kp*error+Ki*integral+Kd*derivative;%output of PID
    d2y = 1 / a * (K * u(k) - b * dy(k) - c * y(k));%y'(t)
    dy(k + 1) = dy(k) + d2y * dt;
    y(k + 1) = y(k) + dy(k + 1) * dt;
    previous_error=error;
end

figure(1)
subplot(3,1,1)
plot(t(1:end-1),e,'linewidth',1.5)
grid
xlabel('time t')
legend({'e(t)'}, 'fontsize', 12)
title(['J = ',num2str(J), ', L = ', num2str(L), ', R = ', num2str(R), ', b = ', num2str(b_fi), ', K = ', num2str(K), ', K_p = ', num2str(Kp), ', K_i = ', num2str(Ki), ', K_d = ', num2str(Kd)])

subplot(3,1,2)
plot(t(1:end-1),u,'linewidth',1.5)
grid
legend({'u(t)'}, 'fontsize', 12)
title('u(t) = ouput of PID ')

subplot(3,1,3)
plot(t,y,'linewidth',1.5)
grid
xlabel('time t')
legend({'\omega(t)'}, 'fontsize', 12)
title('y(t) = output of the system')
