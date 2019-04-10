% ************************************************************************
% * AUTHOR(S) :
% *     Bruno González Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW01.m
% *
% * DESCRIPTION :
% *     Computaciónn Aplicada (Ene 19 Gpo 1)
% *     Homework on Optimization
% *
% * NOTES :
% *     
% *
% * START DATE :
% *     10 Apr 2019
% ************************************************************************

warning('off')
clc;
clear all;
close all;

%% ************************************************************************
% Problem 1:
% Solve the following problem using the optimization toolbox:
%
%            /            2      4 \
%            |     2.1 x_1    x_1  |    2                         2     2
% min f(x) = | 4 - ------- +  ---- | x_1  + x_1 x_2 + (- 4 + 4 x_2 ) x_2
%  x         \         10       3  /
%
% Use function fmincon to solve the problem
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and required results.

objective = (4 - 2.1 * x1^2 + x1^4/3) * x1^2 + x1 * x2 + (-4 + 4 * x2^2) * x2^2;
x0 = [-2,5 2.5]';
disp(['Initial Objective: ' num2str(objective(x0))])

function [c, ceq] = nlcon(x)
    c = 2;
    ceq = 4;
end

A = [];
b = [];
Aeq = [];
beq = [];
lb = [-3 -2]';
ub = [3 2]';
nonlcon = @nlcon;

[x, fval, ef, output, lambda] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon);

disp(x)
disp(['Final Objective: ' num2str(objective(x))])

%% ************************************************************************
% Problem 2:
% Using function fminsearch minimize Branin’s function:
%
%                      2             2
% f(x) = a (x_2 - b x_1  + c x_1 - r)  + s  (1 - t) cos(x_1) + s
%
% where
a = 1;
b = 5.1/(4*pi^2);
c = 5/pi;
r = 6;
s = 10;
t = 1/(8*pi);
%
% for
% x_1 >= -5; x_1 <= 10
% x_2 >=  0; x_2 <= 15
%
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and required results.

% function y = branin(xx)
% % Branin's function
% x1 = xx(1);
% x2 = xx(2);
% t = 1 / (8*pi);
% s = 10;
% r = 6;
% c = 5/pi;
% b = 5.1 / (4*pi^2);
% a = 1;
% term1 = a * (x2 - b*x1.^2 + c*x1 - r).^2;
% term2 = s*(1-t)*cos(x1);
% y = term1 + term2 + s;
% end

% function y = branin2(x1,x2)
% Branin's function (used for plotting)
% t = 1 / (8*pi);
% s = 10;
% r = 6;
% c = 5/pi;
% b = 5.1 / (4*pi^2);
% a = 1;
% term1 = a * (x2 - b*x1.^2 + c*x1 - r).^2;
% term2 = s*(1-t)*cos(x1);
% y = term1 + term2 + s;
% end

%% Solution using fminsearch
% x0 = [rand*15-5 rand*15];
% xmin = fmincon(@branin,x0,[],[],[],[],[-5 0]',[10 15]');
% branin(xmin)

%% Branin's function graph
% [X1,X2] = meshgrid(linspace(-5,10),linspace(0,15));
% F = branin2(X1,X2);
% surf(X1,X2,F)
% shading interp
% xlabel('x_1')
% ylabel('x_2')
