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
[time_2, Y_2] = ode45(@reactub_2, Tspan, Y0, OPTIONS, p);
[time_3, Y_3] = ode45(@reactub_3, Tspan, Y0, OPTIONS, p);
[time_4, Y_4] = ode45(@reactub_4, Tspan, Y0, OPTIONS, p);

% group all data / prepare to plot ...
time     = {time_2, time_3, time_4};
Y        = {   Y_2,    Y_3,    Y_4};
Yprime   = {  Y_2',   Y_3',   Y_4'};
plotName = { 'Oh2',  'Oh3',  'Oh4'};

%% Plot
% Plot limits
noOf_curvesToPlot = 10;
dlim              = 0.02;
time_lim          = [0 - dlim, 1 + dlim];
Y_lim             = [0 - dlim, 1 + dlim];
xi_lim            = [0 - dlim, 1 + dlim];
Yprime_lim        = [0 - dlim, 1 + dlim];

for plotCount = 1:1:3
    % Display concentration vs. time
    totalNoOf_curves  = size(Y{1, plotCount}, 2);
    noOf_curvesToSkip = fix(totalNoOf_curves/noOf_curvesToPlot);
    figure;
    subplot(1, 2, 1)
    for n = linspace(1, totalNoOf_curves, totalNoOf_curves)
        hold all
        if mod(n, noOf_curvesToSkip) == 0
            plot(time{1, plotCount}, Y{1, plotCount}(:, n));
        end
    end
    xlabel('time \tau');
    ylabel('Concentration mol/dm^3');
    axis([time_lim(1) time_lim(2) Y_lim(1) Y_lim(2)])

    % Display concentration vs. distance
    totalNoOf_curves  = size(Yprime{1, plotCount}, 2);
    noOf_curvesToSkip = fix(totalNoOf_curves/noOf_curvesToPlot);
    %figure;
    subplot(1, 2, 2)
    for n = linspace(1, totalNoOf_curves, totalNoOf_curves)
        hold all
        if mod(n, noOf_curvesToSkip) == 0
            plot(xi, Yprime{1, plotCount}(:, n));
        end
    end
    xlabel('distance x/L');
    ylabel('Concentration mol/dm^3');
    axis([xi_lim(1) xi_lim(2) Yprime_lim(1) Yprime_lim(2)])

    % Export Graphics
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 12 6];
    print(plotName{plotCount}, '-dpng', '-r0')
end

%% Print error between O(h)s
sprintf(                                                          ...
    'O(h^2) - O(h^3) error: %f%%',                                ...
    (Yprime{1, 1}(end) - Yprime{1, 2}(end))/Yprime{1, 1}(end)*100 ...
)
sprintf(                                                          ...
    'O(h^2) - O(h^4) error: %f%%',                                ...
    (Yprime{1, 1}(end) - Yprime{1, 3}(end))/Yprime{1, 1}(end)*100 ...
)
sprintf(                                                          ...
    'O(h^3) - O(h^4) error: %f%%',                                ...
    (Yprime{1, 2}(end) - Yprime{1, 3}(end))/Yprime{1, 2}(end)*100 ...
)
