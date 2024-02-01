clear
close all

% GIVEN
m = 1
v_0 = 5
alpha = 65 / 180 * pi
b = 1
delta_t = 0.0001
g = 9.81
%SOLVE

%1) air resistance F(t)=0

dt = 0.01 %time step
s1x(1) = 0 %s(0)
s1y(1) = 0 %s(0)
v1x(1) = v_0 * cos(alpha) %v(0)
v1y(1) = v_0 * sin(alpha) %v(0)
k = 1
while s1y(k) >= 0
      F1x = 0 %air resistance
      F1y = 0
      a1x(k) = 1 / m * F1x;%acceleration
      a1y(k) = 1 / m * (F1y - m * g)

      v1x(k + 1) = v1x(k) + a1x(k) * dt;
      v1y(k + 1) = v1y(k) + a1y(k) * dt;

      s1x(k + 1) = s1x(k) + v1x(k) * dt + 1/2 * a1x(k) * dt ^ 2;
      s1y(k + 1) = s1y(k) + v1y(k) * dt + 1/2 * a1y(k) * dt ^ 2;
      k = k+1;
end

t1=0:dt:(k-1)*dt;
T1=t1(end)
s1=s1x(end)

%2)  F(t)=b*v

dt = 0.01 %time step
s2x(1) = 0 %s(0)
s2y(1) = 0 %s(0)
v2x(1) = v_0 * cos(alpha) %v(0)
v2y(1) = v_0 * sin(alpha) %v(0)
k = 1
while s2y(k) >= 0
      F2x = -b * v2x(k) %air resistance
      F2y = -b * v2y(k)
      a2x(k) = 1 / m * F2x;%acceleration
      a2y(k) = 1 / m * (F2y - m * g)

      v2x(k + 1) = v2x(k) + a2x(k) * dt;
      v2y(k + 1) = v2y(k) + a2y(k) * dt;

      s2x(k + 1) = s2x(k) + v2x(k) * dt + 1/2 * a2x(k) * dt ^ 2;
      s2y(k + 1) = s2y(k) + v2y(k) * dt + 1/2 * a2y(k) * dt ^ 2;
      k = k+1;
end

t2=0:dt:(k-1)*dt;
T2=t2(end)
s2=s2x(end)


%3) F(t)=b*v^2

dt = 0.01 %time step
s3x(1) = 0 %s(0)
s3y(1) = 0 %s(0)
v3x(1) = v_0 * cos(alpha) %v(0)
v3y(1) = v_0 * sin(alpha) %v(0)
k = 1
while s3y(k) >= 0
      F3x = -b * v3x(k) * v3x(k) %air resistance
      F3y = -b * v3y(k) * v3y(k)
      a3x(k) = 1 / m * F3x;%acceleration
      a3y(k) = 1 / m * (F3y - m * g)

      v3x(k + 1) = v3x(k) + a3x(k) * dt;
      v3y(k + 1) = v3y(k) + a3y(k) * dt;

      s3x(k + 1) = s3x(k) + v3x(k) * dt + 1/2 * a3x(k) * dt ^ 2;
      s3y(k + 1) = s3y(k) + v3y(k) * dt + 1/2 * a3y(k) * dt ^ 2;
      k = k + 1;
end

t3=0:dt:(k-1)*dt;
T3=t3(end)
s3=s3x(end)

figure(1)
hold on
plot(s1x,s1y,'r','linewidth',1.5)
plot(s2x,s2y,'green','linewidth',1.5)
plot(s3x,s3y,'blue','linewidth',1.5)

grid
title(['m = ',num2str(m),', v_0 = ',num2str(v_0),', a = ',num2str(alpha), ', b = ',num2str(b), ', alpha = ',num2str(65), ', delta_t = ',num2str(delta_t),
      'T_1 = ',num2str(T1), ', T_2 = ',num2str(T2), ', T_3 = ',num2str(T3), ', s_1 = ',num2str(s1), ', s_2 = ',num2str(s2), ', s_3 = ',num2str(s3),
])


