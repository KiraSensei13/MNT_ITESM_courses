%%
clear;
load('chimneyPDEtoolExports.mat');

chimneyTransientTemperature = u;
chimneyMeshEdges            = e;
chimneyMeshPoints           = p;
chimneyMeshTriangles        = t;
chimneyMovie                = M;

%%
%     +---------------- boundary w ----------------+
%     |            +--- boundary y ---+            |
%     |            |                  |            |
% boundary 2   boundary 4         boundary z   boundary x
%     |            |                  |            |
%     |            +--- boundary 3 ---+            |
%     +---------------- boundary 1 ----------------+
% Due to symetry,
% boundary 1 = boundary w
% boundary 2 = boundary x
% boundary 3 = boundary y
% boundary 4 = boundary z

% Specifying the coordinates of the points at boundary 1
nx1 = 100; % number of increments in boundary 1
xb1 = linspace(0.00, 4.00, nx1 + 1);
yb1 = 0.00;

% Specifying the coordinates of the points at boundary 2
ny2 = 100; % number of increments in boundary 2
xb2 = 0.00;
yb2 = linspace(0.00, 1.00, ny2 + 1);

% Specifying the coordinates of the points at boundary 3
nx3 = 100; % number of increments in boundary 3
xb3 = linspace(0.20, 3.80, nx3 + 1);
yb3 = 0.20;

% Specifying the coordinates of the points at boundary 4
ny4 = 100; % number of increments in boundary 4
xb4 = 0.20;
yb4 = linspace(0.20, 0.80, ny4 + 1);

% Increment of distance between points at each boundary
dx1 = distanceIncrement(xb1, nx1);
dy2 = distanceIncrement(yb2, ny2);
dx3 = distanceIncrement(xb3, nx3);
dy4 = distanceIncrement(yb4, ny4);

% External and internal heat transfer coefficients [W/m2-K]
ho = 10;
hi = 50;

% External and internal temperatures
To = 25;
Ti = 500;

% Loop to calculate heat rate per unit length at each time
[n, dt] = size(u);                % dt, time increment specified in PDEtool
time   = linspace(0, 360000, dt); % time domain in seconds
time   = time / 3600;             % time domain in hours
To_avg = zeros(dt, 1);
Ti_avg = zeros(dt, 1);
Qo     = zeros(dt, 1);
Qi     = zeros(dt, 1);
for i = 1:dt
    T = u(:, i);
    [Tb1, Qb1] = Q_boundery(p, t, T, xb1, yb1,  ho, dx1, To);
    [Tb2, Qb2] = Q_boundery(p, t, T, xb2, yb2,  ho, dy2, To);
    [Tb3, Qb3] = Q_boundery(p, t, T, xb3, yb3, -hi, dx3, Ti);
    [Tb4, Qb4] = Q_boundery(p, t, T, xb4, yb4, -hi, dy4, Ti);
    To_avg(i) = (Tb1 + Tb2) / 2;
    Ti_avg(i) = (Tb3 + Tb4) / 2;
    Qo(i)     = 2 * (Qb1 + Qb2);
    Qi(i)     = 2 * (Qb3 + Qb4);
end

%%
% plot heat rate
figure
plot(time, (abs(Qo)), ...
     time, (abs(Qi)))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('external', 'internal')
xlabel("time [h]")
ylabel("heatRate/L [W/m]")

% plot temperature
figure
plot(time, (abs(To_avg)), ...
     time, (abs(Ti_avg)))
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
legend('external', 'internal')
xlabel("time [h]")
ylabel("temperature [K]")

% animate temperature
figure
movie(gcf, M)

% save animation as .gif
movie2gif( M, 'animation.gif')

%%
function dxy = distanceIncrement(xyb, nxy)
    dxy = (max(xyb) - min(xyb)) / nxy;
end

function [Tavg, Qb] = Q_boundery(p, t, T, xb, yb, h, dxy, Tref)
    Tb   = tri2grid(p, t, T, xb, yb);
    Tavg = mean(Tb);
    Qb   = h * dxy * (sum(Tb - Tref) - (1/2) * (Tb(1) + Tb(end) - 2 * Tref));
end