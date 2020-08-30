%% Setup
clear;
close all;

%% Import data
FEATool_exports = csvread('./FEATool.csv');

%% Plotting
N = size(FEATool_exports, 2)/2; % number of curves

% plot limits
dlim = 0.02;
x_lim = [0 - dlim, 1 + dlim];
t_lim = [0 - dlim, 1 + dlim];
u_lim = [0 - dlim, 1 + dlim];

% 2D line plot
figure
for n = linspace(1, N, N)
    data = FEATool_exports(:, n*2 - 1:n*2);
    data = sortrows(data,1);
    size(data)
    x = data(:, 1);
    u = data(:, 2);
%     x = FEATool_exports(:, n*2 - 1);
%     u = FEATool_exports(:, n*2);
    hold all
    plot(x, u)
end
xlabel('Distance x')
ylabel('Species u')
axis([x_lim(1) x_lim(2) u_lim(1) u_lim(2)])