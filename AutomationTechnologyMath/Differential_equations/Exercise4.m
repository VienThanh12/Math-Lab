clear
close all

M = 100
m = 10
k = 500000
b = 100
R = 0.1
w_given = 100
% GIVEN

% My′′(t)+by′(t)+ky(t)=m R ω^2 sin(ωt), y(0) = 0,y′(0) = 0
% a*y''(t)+b*y'(t)+c*y(t)=A*sin(w*t+theta)+B, y(0)=y0, x'(0)=y1


a = M
%b = b
c = k
A = m * R * w_given * w_given
B = 0
y0h = 0
y1h = 0

K=1/sqrt((b*w_given)^2+(c-a*w_given^2)^2)
phi=atan2(-b*w_given,c-a*w_given^2)

xp=@(t) K*A*sin(w_given*t+phi)+B/c


y0 = y0h - K*A*sin(phi)-B/c
y1 = y1h - K*A*w_given*cos(phi)

D = b ^ 2 - 4 * a * c

if D<0
    alpha=-b/(2*a);
    tau=-1/alpha;
    w=sqrt(-D)/(2*a);
    C1=y0;
    C2=(y1-alpha*y0)/w;
    C=sqrt(C1^2+C2^2);
    th=atan2(C1,C2);
    y=@(t) C*exp(alpha*t).*sin(w*t+th);
    dy=@(t) C*exp(alpha*t)*alpha.*sin(w*t+th)+C*exp(alpha*t).*cos(w*t+th)*w;
    d2y=@(t) C*exp(alpha*t)*alpha^2.*sin(w*t+th) +2*C*exp(alpha*t)*alpha.*cos(w*t+th)*w-C*exp(alpha*t).*sin(w*t+th)*w^2;

elseif D==0
    la=-b/(2*a);
    tau=-1/la;
    C1=y0;
    C2=y1+la*y0;
    y=@(t) C1*exp(la*t)+C2*t.*exp(la*t);
    dy=@(t) C1*exp(la*t)*la+C2*exp(la*t)+C2*t.*exp(la*t)*la;
    d2y=@(t) C1*exp(la*t)*la^2+2*C2*la*exp(la*t) + C2*t.*exp(la*t)*la^2;


else
    la1=(-b+sqrt(D))/(2*a);
    la2=(-b-sqrt(D))/(2*a);
    tau=-1/la1;
    C1=(y1-la2*y0)/(la1-la2);
    C2=(y1-la1*y0)/(la2-la1);
    y=@(t) C1*exp(la1*t)+C2*exp(la2*t);
    dy=@(t) C1*exp(la1*t)*la1+C2*exp(la2*t)*la2;
    d2y=@(t) C1*exp(la1*t)*la1^2+C2*exp(la2*t)*la2^2;
end

x = @(t) xp(t) + y(t)

tmax=tau
t=0:tmax/1000:1;

figure(1)
plot(t, x(t), 'linewidth', 1)
grid
title(['M = ',num2str(M),', m = ',num2str(m),', k = ',num2str(k),...
      ', b = ',num2str(b),', R = ',num2str(R),', \omega   = ',num2str(w_given)])
xlabel('aika t')
ylabel('paikka y(t)')

K_w = @(w) 1 ./ sqrt((b .* w) .^ 2 + (c - a .* (w .^ 2)) .^2)
A_w = @(w) m .* R .* w .* w

w = 0:0.1:10 * sqrt(k / M)
figure(2)
hold
plot(100, A_w(100) .* K_w(100), 'r.', 'markersize', 15)
plot(w, A_w(w) .* K_w(w), 'blue', 'linewidth', 1.5)
grid
xlabel('kulmanopeus \omega')
ylabel('y(t ):n amplitudi')
