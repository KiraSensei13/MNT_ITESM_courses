% Program to se the velocity field and to track a tracer figure
% By JLLopez CFD 09/29/2020

% HW07
% adapted from (jose lopez salinas)'s solution 
clear;
close all;

% The tracer figure is a circle
rho     = 0.3;           % Radius of the circle
x0      = 0.5;
y0      = 0.5;           % Center of the circle
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

YL1 = YS(round(size(YS,1)/3),:); % get the time position at 1/3 of the time
YL2 = YS(round(size(YS,1)*2/3),:); % get the time position at 2/3 of the time
YLF = YS(end,:); % get the last time position
for i = 1 : m
    x1(i) = YL1(2 * i - 1);
    y1(i) = YL1(2 * i);
    x2(i) = YL2(2 * i - 1);
    y2(i) = YL2(2 * i);
    xf(i) = YLF(2 * i - 1);
    yf(i) = YLF(2 * i);
end

% plots arrows with directional components U and V at the Cartesian coordinates specified by X and Y
figure;
hold all;
[xx, yy, Ux, Uy] = MPlotxx1();

% plots initial position and final to compare
plot( ...
    x,   y, 'ro-', ...
    x1, y1, 'go-', ...
    x2, y2, 'bo-', ...
    xf, yf, 'mo-');
xlabel('x');
ylabel('y');
%xlim([0.2 1.4]);ylim([0.2 1.4]);

% Export Graphics
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 6];
print('velocityField', '-dpng', '-r0')