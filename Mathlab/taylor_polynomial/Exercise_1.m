clear
close all
#GIVEN

n = 5

#SOLVE

Tn = 0
x = -1:0.01:1
for k = 1:n
   Tn = Tn + (-1)  .^ (k - 1) ./ k .* x .^ (k)
end

#graph
hold on
f = @(x) log(x + 1)
plot(x, f(x), 'blue' ,'linewidth', 1.5)
plot(x,Tn,'r','linewidth',1.2)
grid
ylim([-2, 2])

title(['n = ', num2str(n)])
xlabel('x')

hold off
