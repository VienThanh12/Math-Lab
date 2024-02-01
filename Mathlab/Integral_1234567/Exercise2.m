clear
close all

% GIVEN

A = 5
T = 4

T = T / 2
% Cal
% u(t) = A / T * x

% Integrate THE MEAN
% The mean
uavg = 1 / T * A .* T ./ 2


% Integrate THE ROOT MEAN SQUARE
% THE ROOT MEAN SQUARE
urms = sqrt(1 / (T) * (A .^ 2) .* T / 3)

