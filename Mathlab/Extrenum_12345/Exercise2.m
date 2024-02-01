clear
close all

v = 50
e = 100
r = 10000

c = @(q) v .* (q ./ 2) + e .* (r ./ q)

% dc = v/2-(e*r)/q^2
% q = sqrt(2)*sqrt(e)*sqrt(r/v) or -sqrt(2)*sqrt(e)*sqrt(r/v)
% q > 0

u = sqrt(2)*sqrt(e)*sqrt(r/v)
c(u)
hold on
q = 0:10:800

%[cmin,indmin] = min(c(q)) %min value of vector f is its element number indmin
%qmin = q(indmin)

plot(q, c(q), 'linewidth', 1.5)
plot(u,c(u),'r.','markersize',15);
xlabel('eräkoko q')
title(['r = ', num2str(r), ', kv = ', num2str(v), ', ke = ', num2str(e), ': min = ', num2str(c(u)), ', kun q = ', num2str(u)])
ylabel('vuotuiset kustammukset')
grid
hold off
