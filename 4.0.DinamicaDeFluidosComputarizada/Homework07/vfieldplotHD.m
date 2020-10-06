% Program to se the velocity field and to track a tracer figure
% By JLLopez CFD 09/29/2020

% HW07
% adapted from (jose lopez salinas)'s solution 
clear;
close all;

% The tracer figure is a circle
rho     = 0.2;         % Radius of the circle
x0      = 0.5;
y0      = 0.75;          % Center of the circle
p1      = [x0 y0 rho];	 % parameter to draw the circle
[x,y]   = ftriangle(p1); % function to generate the shape
[n1,m1] = size(x);
m       = max(n1,m1);
to      = 0;
tf      = 1.00;

% vector position
for i = 1 : m
    z(2 * i - 1) = x(i);
    z(2 * i)     = y(i);
end

% parameters you may need in the vector field function
p     = 1;
zo    = z;                  % initial condition
tspan = linspace(0,0.5,20);	% time span to track the fluid parcels
tf    = [0 1];

% Solution of the ODEs dr/dt , here you solve the velocity field eqn (vector field) 
[time, YS] = ode45(@Vfield, tspan, zo, [], p);

% You will plot that last time position
YL=YS(end,:);
for i = 1 : m
    xf(i) = YL(2*i-1);
    yf(i) = YL(2*i);
end

% p3      = p1;
% [xs,ys] = fsquarex(p3);
% plot(xs,ys,'rs-');
% hold;

% plots arrows with directional components U and V at the Cartesian coordinates specified by X and Y
hold;
[x2, y2, Ux, Uy] = MPlotxx1();

% plots initial position and final to compare
plot(x, y, 'bo-', xf, yf, 'ro-');
xlabel('x');
ylabel('y');
%xlim([0.2 1.4]);ylim([0.2 1.4]);