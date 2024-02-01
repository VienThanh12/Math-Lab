clear
close all

x=@(t) 1.5*sin(0.5*t)
y=@(t) cos(1.25*t)
vx=@(t) 0.75*cos(0.5*t)
vy=@(t) -1.25*sin(1.25*t)
v=@(t) sqrt(vx(t).^2+vy(t).^2)
ax=@(t) -0.375*sin(0.5*t)
ay=@(t) -1.5625*cos(1.25*t)
a= @(t) sqrt(ax(t).^2+ay(t).^2)
aT= @(t) (vx(t).*ax(t)+vy(t).*ay(t))./v(t)
aN= @(t) (vx(t).*ay(t)-vy(t).*ax(t))./v(t)



T=4*pi
t=0:T/1000:T;
t0=7.4



%% animation

plot(x(t),y(t),'linewidth',1.5)



