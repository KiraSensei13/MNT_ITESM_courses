%% HW05 part B
%% TUBULAR REACTOR adapted from (jose lopez salinas)'s solution 
% Setup
clear;
close all;

%% Runge-Kutta
% Use any integration technique (Runge-Kutta type scheme) to solve in time,
% and finite differences to integrate in axial domain.

p(1)  = 0.001; % Diffusion coefficient D
p(2)  = 1.0;    % Injection concentration c0
p(3)  = 1.5;    % First order kinetic coefficient k
p(4)  = 1.0;    % Velocity of fluid injection vo
M     = 2*640;  % Number of nodes
p(5)  = M;
Tspan = [0 1];  % Domain of time
xi    = linspace(0, 1, M);

% Initial conditions of the resulting set of ODEs
Y0    = zeros(M, 1);
Y0(1) = 1.0;

% Solve differential equation (medium order method)
% use @reactub_2 for O(h^2) truncation error
% use @reactub_3 for O(h^3) truncation error
% use @reactub_4 for O(h^4) truncation error
OPTIONS   = [];
[time, Y] = ode45(@reactub_2, Tspan, Y0, OPTIONS, p);
Yprime    = Y';

% plot limits
noOf_curvesToPlot = 25;
dlim              = 0.02;
time_lim          = [min(time)        - dlim, max(time)        + dlim];
Y_lim             = [min(min(Y))      - dlim, max(max(Y))      + dlim];
xi_lim            = [min(xi)          - dlim, max(xi)          + dlim];
Yprime_lim        = [min(min(Yprime)) - dlim, max(max(Yprime)) + dlim];

% Display
totalNoOf_curves  = size(Y, 2);
noOf_curvesToSkip = fix(totalNoOf_curves/noOf_curvesToPlot);
figure;
for n = linspace(1, totalNoOf_curves, totalNoOf_curves)
    hold all
    if mod(n, noOf_curvesToSkip) == 0
        plot(time, Y(:, n))
    end
end
xlabel('\tau');
ylabel('Concentration mol/dm^3');
axis([time_lim(1) time_lim(2) Y_lim(1) Y_lim(2)])

% Display
totalNoOf_curves  = size(Yprime, 2);
noOf_curvesToSkip = fix(totalNoOf_curves/noOf_curvesToPlot);
figure;
for n = linspace(1, totalNoOf_curves, totalNoOf_curves)
    hold all
    if mod(n, noOf_curvesToSkip) == 0
        plot(xi, Yprime(:, n))
    end
end
xlabel('distance x/L');
ylabel('Concentration mol/dm^3');
axis([xi_lim(1) xi_lim(2) Yprime_lim(1) Yprime_lim(2)])