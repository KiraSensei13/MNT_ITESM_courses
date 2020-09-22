%% HW05 part C
%% GROWING BUBBLES
% Setup
clear;
close all;

% Runge-Kutta
% Use any integration technique (Runge-Kutta type scheme) to solve in time,
% and finite differences to integrate in axial domain.

% k = 1.4 adiabatic process, k = 1 isothermic
% alphaM = 0 inviscid
% betaM = 0 negligible surface tension
k      = {1.0, 1.4, 1.7};
alphaM = {0.0, 0.5, 0.1};
betaM  = {0.0, 0.2, 0.7};

% Plot
figure;
solveNplot_growingBubbles(k{1}, alphaM{1}, betaM{1}, sprintf('k = %1.1f', k{1}))
solveNplot_growingBubbles(k{2}, alphaM{1}, betaM{1}, sprintf('k = %1.1f', k{2}))
solveNplot_growingBubbles(k{3}, alphaM{1}, betaM{1}, sprintf('k = %1.1f', k{3}))
% Export Graphics
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 8];
suptitle(sprintf('Gas Molecule Shape and Size Effect : \\alpha = %1.1f and \\beta = %1.1f', alphaM{1}, betaM{1}))
print('growingBubbles_GasMoleculeShapeAndSizeEffect_Alpha0Beta0', '-dpng', '-r0')

figure;
solveNplot_growingBubbles(k{2}, alphaM{1}, betaM{1}, sprintf('\\alpha = %1.1f', alphaM{1}))
solveNplot_growingBubbles(k{2}, alphaM{2}, betaM{1}, sprintf('\\alpha = %1.1f', alphaM{2}))
solveNplot_growingBubbles(k{2}, alphaM{3}, betaM{1}, sprintf('\\alpha = %1.1f', alphaM{3}))
% Export Graphics
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 8];
suptitle(sprintf('Effect of the Viscosity : k = %1.1f and \\beta = %1.1f', k{2}, betaM{1}))
print('growingBubbles_EffectOfTheViscosity_K1_4Beta0', '-dpng', '-r0')

figure;
solveNplot_growingBubbles(k{2}, alphaM{1}, betaM{1}, sprintf('\\beta = %1.1f', betaM{1}))
solveNplot_growingBubbles(k{2}, alphaM{1}, betaM{2}, sprintf('\\beta = %1.1f', betaM{2}))
solveNplot_growingBubbles(k{2}, alphaM{1}, betaM{3}, sprintf('\\beta = %1.1f', betaM{3}))
% Export Graphics
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 8];
suptitle(sprintf('Effect of the Surface Tension : k = %1.1f and \\alpha = %1.1f', k{2}, alphaM{1}))
print('growingBubbles_EffectOfTheViscosity_K1_4Alpha0', '-dpng', '-r0')

%% DRAINING TANK
% Setup
clear;
close all;

% Runge-Kutta
% Use any integration technique (Runge-Kutta type scheme) to solve in time,
% and finite differences to integrate in axial domain.



function[] = solveNplot_growingBubbles(k, alphaM, betaM, curveLabel)
    tspan = linspace(0, 35, 500);
    y1    = 1;
    y2    = 0;
    yo    = [y1, y2]';

    % Solve differential equation (medium order method)
    par(1) = k;
    par(2) = alphaM;
    par(3) = betaM;
    [t,Y]  = ode45(@growingBubbles, tspan, yo, [], par);
    time   = t;
    Yout   = Y;

    %% Plot
    % Display radius vs. time
    subplot(2,1,1);
    hold all
    plot(time, Yout(:,1), 'DisplayName', curveLabel);
    xlabel('\tau');ylabel('R/Ro');
    legend

    % Display d(radius) vs. time
    subplot(2,1,2);
    hold all
    plot(time, Yout(:,2), 'DisplayName', curveLabel);
    xlabel('\tau');ylabel('d(R/Ro) /d\tau');
    legend
end