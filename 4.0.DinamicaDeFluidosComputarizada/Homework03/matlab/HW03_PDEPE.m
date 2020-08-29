%% Setup
clear;
close all;

%% Runge-Kutta
% Use any integration technique (Runge-Kutta type scheme) to solve in time,
% and finite differences to integrate in axial domain.

p(1)  = 0.01; % Diffusion coefficient D
p(2)  = 1.0;  % Injection concentration c0
p(3)  = 0.5;  % First order kinetic coefficient k
p(4)  = 1.0;  % Velocity of fluid injection vo
M     = 10;   % Number of nodes
p(5)  = M;
Tspan = [0 1]; % Domain of time

% Initial conditions of the resulting set of ODEs
Y0    = zeros(M,1);
Y0(1) = 1.0;

% Solve differential equation (medium order method)
OPTIONS = [];
[time, Y] = ode45(@reactub, Tspan, Y0, OPTIONS, p);

% Display
figure
plot(time, Y)
xlabel('Time t')
ylabel('Injection concentration C')

%% Diffusion
% This is the main function. Within this function the meshes are defined,
% PDEPE is called and the results are plotted

% Parameters
% dC/dt=div(D Grad C)-k*C- dot(v,Grad C)
% dC/dt= D d(dC/dx)/dx- k*C-v*dC/dx
% in the previous pde, all are partial derivatives
P(1) = 0.01;                   % Diffusion coefficient D
P(2) = 1.0;                    % Injection concentration c0
P(3) = 0.5;                    % First order kinetic coefficient k
P(4) = 1.0;                    % Velocity of fluid injection vo
L    = 1;                      % Length of domain
maxt = 1;                      % Max. simulation time
m    = 0;                      % Parameter corresponding to the symmetry of
%                                the problem (see help)
t    = linspace(0, maxt, 100); % Tspan
x    = linspace(0, L, 100);    % xmesh

% DiffusionPDEfun: Function containg the PDEs
% DiffusionICfun: Function containing the ICs for t = 0 at all x
% DiffusionBCfun: Function containing the BCs for x = 0 and x = L
% PDEPE returns the solution as multidimensional array of size
% xmesh x tspan x (# of variables)
sol = pdepe(          ...
    m,                ...
    @DiffusionPDEfun, ... % Function containg the PDEs
    @DiffusionICfun,  ... % Function containing the ICs for t = 0 at all x
    @DiffusionBCfun,  ... % Function containing the BCs for x = 0 and x = L
    x,                ...
    t,                ...
    [],               ... % Options
    P                 ... % Parameters
);
u = sol;

% Plotting
% 3D surface plot
figure
surf(x, t, u, 'edgecolor', 'none');
xlabel('Distance x')
ylabel('Time t')
zlabel('Species u')
axis([0 L 0 maxt 0 P(2)])
set(gcf(), 'Renderer', 'painters')

% 2D line plot
figure
hold all
for n = linspace(1, length(t), 10)
    plot(x, sol(n, :))
end
xlabel('Distance x')
ylabel('Species u')
axis([0 L 0 P(2)])