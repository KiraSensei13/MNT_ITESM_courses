%% TUBULAR REACTOR adapted from (jose lopez salinas)'s solution 
% Setup
clear;
close all;

%% Runge-Kutta
% Use any integration technique (Runge-Kutta type scheme) to solve in time,
% and finite differences to integrate in axial domain.

p(1)  = 0.01;  % Diffusion coefficient D
p(2)  = 1.0;   % Injection concentration c0
p(3)  = 0.5;   % First order kinetic coefficient k
p(4)  = 1.0;   % Velocity of fluid injection vo
M     = 10;    % Number of nodes
p(5)  = M;
Tspan = [0 1]; % Domain of time

% Initial conditions of the resulting set of ODEs
Y0    = zeros(M,1);
Y0(1) = 1.0;

% Solve differential equation (medium order method)
OPTIONS = [];
[time, Y] = ode45(@reactub, Tspan, Y0, OPTIONS, p);

% Display
% figure
% plot(time, Y)
% xlabel('Time t')
% ylabel('Injection concentration C')

%% PDEPE - Diffusion
% This is the main function. Within this function the meshes are defined,
% PDEPE is called and the results are plotted

% Parameters
% dC/dt = div(D Grad C) - k*C - dot(v, Grad C)
% dC/dt = D d(dC/dx)/dx - k*C - v*dC/dx
% in the previous pde, all are partial derivatives
P(1) = 0.01;                    % Diffusion coefficient D
P(2) = 1.0;                     % Injection concentration c0
P(3) = 0.5;                     % First order kinetic coefficient k
P(4) = 1.0;                     % Velocity of fluid injection vo
L    = 1;                       % Length of domain
maxt = 1;                       % Max. simulation time
m    = 0;                       % Parameter corresponding to the symmetry
%                                 of the problem (0 for slab, 1 for
%                                 cylinder, or 2 for sphere)
step = 32;
t    = linspace(0, maxt, step); % Tspan
x    = linspace(0, L, step);    % xmesh

% PDEPE returns the solution as multidimensional array of size
% xmesh x tspan x (# of variables)
sol = pdepe(          ...
    m,                ...
    @DiffusionPDEfun, ... % Function containing the PDEs
    @DiffusionICfun,  ... % Function containing the ICs for t=0 at all x
    @DiffusionBCfun,  ... % Function containing the BCs for x=0 and x=L
    x,                ... % Spatial mesh
    t,                ... % Time span of integration
    [],               ... % Options
    P                 ... % Parameters
);

% The element ui(j,k) = sol(j,k,i) approximates ui at
% ui(t,x) = (tspan(j),xmesh(k)).
u = sol;

% Plotting
% plot limits
dlim = 0.02;
x_lim = [0 - dlim, L + dlim];
t_lim = [0 - dlim, maxt + dlim];
u_lim = [0 - dlim, P(2) + dlim];

% 3D surface plot
figure
surf(x, t, u, 'edgecolor', 'none');
xlabel('Distance x')
ylabel('Time t')
zlabel('Species u')
axis([x_lim(1) x_lim(2) t_lim(1) t_lim(2) u_lim(1) u_lim(2)])
set(gcf(), 'Renderer', 'painters')

% 2D line plot
figure
for n = linspace(1, step, step)
    hold all
    plot(x, u(n, :))
end
xlabel('Distance x')
ylabel('Species u')
axis([x_lim(1) x_lim(2) u_lim(1) u_lim(2)])

% 2D line plot
figure
for n = linspace(1, step, step)
    hold all
    plot(t, u(:, n))
end
xlabel('Time t')
ylabel('Species u')
axis([t_lim(1) t_lim(2) u_lim(1) u_lim(2)])

%% GUI script steps of HW03_FEATool.fes
% 
% % INITIALIZATION AND MODEL SELECTION
% 1. To start a new model click the New Model toolbar button, or select
% New Model... from the File menu.
% 2. Select the 1D radio button.
% 3. Select the Convection and Diffusion physics mode from the Select
% Physics drop-down menu.
% 4. Press OK to finish the physics mode selection.
% 
% % GEOMETRY
% 5. Press the Create line Toolbar button.
% 6. Press OK to finish and close the dialog box.
% 
% % DRID
% 7. Switch to Grid mode by clicking on the corresponding Mode Toolbar
% button.
% 8. Enter  0.01  into the Grid Size edit field.
% 
% % EQUATION
% 9. Switch to Equation mode by clicking on the corresponding Mode
% Toolbar button.
% 10. Enter  0.1  into the Time scaling coefficient edit field.
% 11. Enter  0.01  into the Diffusion coefficient edit field.
% 12. Enter  1  into the Convection velocity in x-direction edit field.
% 13. Enter  0.0  into the Reaction rate edit field.
% 14. Press the Apply button.
% 15. Press OK to finish the equation and subdomain settings
% specification.
% 
% % BOUNDARY CONDITIONS
% 16. Switch to Boundary mode by clicking on the corresponding Mode
% Toolbar button.
% 17. Select Concentration from the Convection and Diffusion drop-down
% menu.
% 18. Enter  1  into the Concentration edit field.
% 19. Select 2 in the Boundaries list box.
% 20. Press the Apply button.
% 21. Press OK to finish the boundary condition specification.
% 
% % SOLVE SETTINGS
% 22. Switch to Solve mode by clicking on the corresponding Mode
% Toolbar button.
% 23. Press the Settings Toolbar button.
% 24. Select Time-Dependent from the Solution and solver type
% drop-down menu.
% 25. Enter  30  into the Maximum number of non-linear iterations edit
% field.
% 26. Enter  0.001  into the Time step size edit field.
% 27. Enter  0.001  into the Duration of time-dependent simulation (maximum
% time) edit field.
% 28. Press the Solve button.
% 
% % SOLVE FOR DIFFERENT TIMES
% 29. Switch to Solve mode by clicking on the corresponding Mode
% Toolbar button.
% 30. Press the Settings Toolbar button.
% 31. Enter  0.010  into the Duration of time-dependent simulation (maximum
% time) edit field.
% 32. Press the Solve button.
% 
% % repeat steps 29. through 32 for different values of "Duration of
% % time-dependent simulation (maximum time)" from 0.020 to 0.140 with 0.010
% % increments
% 
% % GUI script HW03_FEATool done!

% Import data
FEATool_exports = csvread('./FEATool.csv');

% Plotting
N = size(FEATool_exports, 2)/2; % number of curves

% plot limits
dlim = 0.02;
x_lim = [0 - dlim, 1 + dlim];
t_lim = [0 - dlim, 1 + dlim];
u_lim = [0 - dlim, 1 + dlim];

% 2D line plot
% for time t = [0.001, 0.140]
figure
for n = linspace(1, N, N)
    data = FEATool_exports(:, n*2 - 1:n*2);
    data = sortrows(data,1);
    x = data(:, 1); % distance
    u = data(:, 2); % species/concentration
    hold all
    plot(x, u)
end
xlabel('Distance x')
ylabel('Concentration C')
axis([x_lim(1) x_lim(2) u_lim(1) u_lim(2)])