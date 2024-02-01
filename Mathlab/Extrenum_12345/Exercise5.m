clear
close all

m = 2
R = 1
k = 3
M = 4
b = 5

A = @(x) (m .* R .* (x .^ 2)) ./ sqrt(((k - M .* (x .^ 2)) .^ 2) + ((b .* x) .^ 2))
w = (sqrt(2) * k) / sqrt(2 * M * k - b ^ 2)
A(w)

x = 0:0.1:10
f = (m .* R .* (x .^ 2)) ./ sqrt(((k - M .* (x .^ 2)) .^ 2) + ((b .* x) .^ 2))


figure(1)
hold on


%[fmax,indmax] = max(f) %largest value of vector f is its element number indmin
%xmax=x(indmax)
plot(x, f, 'linewidth',1.5)
ylabel('A(w)')
xlabel('w')

if(2 * M * k - b ^ 2 < 0)
  title(['m = ', num2str(m), ', R = ', num2str(R), ', k = ', num2str(k), ', M = ', num2str(M), ', b = ', num2str(b)])
else {
  plot(w, A(w), 'r.','markersize', 10)
  title(['m = ', num2str(m), ', R = ', num2str(R), ', k = ', num2str(k), ', M = ', num2str(M), ', b = ', num2str(b), ' : max A(', num2str(w), ')', ' = ', num2str(A(w))])
  }
end
grid
hold off
