clear
close all

%GIVEN
T1 = 1
T2 = 2
T3 = 3
amax = 2
T = 4 * T1 + 2 * T2 + T3

t1 = T1
t2 = T1 + T2
t3 = T1 + T2 + T1
t4 = T1 + T2 + T1 + T3
t5 = T1 + T2 + T1 + T3 + T1
t6 = T1 + T2 + T1 + T3 + T1 + T2
t7 = T1 + T2 + T1 + T3 + T1 + T2 + T1


dt = T / 1000
t = 0:dt:T;


a1 = amax * t
a2 = amax
a3 = -amax * t + 8
a4 = 0
a5 = -amax * t + 14
a6 = -amax
a7 = amax * t - 22

N=length(t)
vk=zeros(1,N);
sk=zeros(1,N);
vk(1)=0;
sk(1)=0;

a = a1.*(t<=t1) + a2.*((t>t1)&(t<t2)) + a3.*((t>t2)&(t<t3)) + a4.*((t>t3)&(t<t4)) + a5.*((t>t4)&(t<t5)) + a6.*((t>t5)&(t<t6)) + a7.*((t>t6)&(t<t7));

for n=1:N-1
    vk(n+1)=vk(n)+a(n)*dt;
    sk(n+1)=sk(n)+vk(n)*dt+1/2*a(n)*dt^2;
end

figure(4)
subplot(3,1,1)
plot(t,a,'linewidth',1.5)
xlim([0, T])
grid
title('kiihtyvyys a(t)')
subplot(3,1,2)
plot(t,vk,'linewidth',1.5)
xlim([0, T])
title('nopeus v(t)')
grid
subplot(3,1,3)
plot(t,sk,'linewidth',1.5)
xlim([0, T])
grid
title('paikka s(t)')
xlabel('aika t')

