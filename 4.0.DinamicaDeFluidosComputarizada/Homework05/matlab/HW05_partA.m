%% HW05 part A - Velocity Field
% adapted from (jose lopez salinas)'s solution 
% Setup
clear;
close all;

% create points to visualize
xyLim  = 2.5;
xyStep = xyLim/10;
[x, y] = meshgrid(-xyLim : xyStep : xyLim);

VectorX = cos(y); % vector in the x direction
VectorY = sin(x); % vector in the y direction

V        = sqrt(VectorX.^2 + VectorY.^2);
PHI      = 6 + x.^3 / 3 - y.^2 .* x - y;
[Dx, Dy] = gradient(V, 0.2, 0.2);

% Display
figure;
quiver(x, y, VectorX, VectorY);
hold on;
contour(x, y, PHI);
colorbar;
hold off;
xlabel('x-axis');
ylabel('y-axis');
title('Velocity Field, and Pressure Lines');

% Export Graphics
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 9 6];
print('vectorField', '-dpng', '-r0')