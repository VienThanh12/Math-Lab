clear
close all

%GIVEN

R = 1
L = 2
C_1 = 0.5
U = 1

% Solution

% L·q′′(t)+R·q′(t)+1/C·q(t)=U
% a*y''(t)+b*y'(t)+c*y(t)=A, y(0)=y0, x'(0)=y1
% solution_formulas.pdf, 2 b i), p.5

a = L
b = R
c = 1 / C_1
A = U

y0h = 0
y1h = 0

y0 = y0h - A / c
y1 = y1h

D = b ^ 2 - 4 * a * c

if D < 0
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

elseif D == 0
    la=-b/(2*a);
    tau=-1/la;
    C1=y0;
    C2=y1+la*y0;
    y=@(t) C1*exp(la*t)+C2*t.*exp(la*t);
    dy=@(t) C1*exp(la*t)*la+C2*exp(la*t)+C2*t.*exp(la*t)*la;
    d2y=@(t) C1*exp(la*t)*la^2+2*C2*la*exp(la*t) + C2*t.*exp(la*t)*la^2;


elseif D > 0
    la1=(-b+sqrt(D))/(2*a);
    la2=(-b-sqrt(D))/(2*a);
    tau=-1/la1;
    C1=(y1-la2*y0)/(la1-la2);
    C2=(y1-la1*y0)/(la2-la1);
    y=@(t) C1*exp(la1*t)+C2*exp(la2*t);
    dy=@(t) C1*exp(la1*t)*la1+C2*exp(la2*t)*la2;
    d2y=@(t) C1*exp(la1*t)*la1^2+C2*exp(la2*t)*la2^2;
end

w = @(t) A / c + y(t)

tmax = 10*tau
t = 0:tmax/1000:tmax;

figure(1)
hold
plot(t, R * (dy(t)),'red','linewidth',1.5)
plot(t, L * (d2y(t)),'green','linewidth',1.5)
plot(t, w(t) * c,'blue','linewidth',1.5)
grid
xlabel('time t')
legend({'U_R','U_L', 'U_C'},'fontsize', 12)
title(['R = ', num2str(R),', L = ', num2str(L), ', C = ', num2str(c), ', U = ', num2str(U)])
