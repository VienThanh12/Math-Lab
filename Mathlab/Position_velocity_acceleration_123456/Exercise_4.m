clear
close all

% GIVEN


r = 2
V = 3
T = 4.1888
t = 0.3 * T


% SOLVE

w = V * t / r

x = @(t) r * cos(w * t)
y = @(t) r * sin(w * t)
T=4*pi
t=0:T/1000:2 * pi;
t0=7.4



%% animation
figure(2)
%for t0=0:T/1000:T
t
plot(x(t), y(t),'linewidth',1.5)
%hold


%hold off
%grid
%axis([-2,2,-2,2])
%axis square
%pause(0.01)
%end
