%% Setup
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

%% Diffusion
% This is the main function. Within this function the meshes are defined,
% PDEPE is called and the results are plotted

% Parameters
% dC/dt=div(D Grad C)-k*C- dot(v,Grad C)
% dC/dt= D d(dC/dx)/dx- k*C-v*dC/dx
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