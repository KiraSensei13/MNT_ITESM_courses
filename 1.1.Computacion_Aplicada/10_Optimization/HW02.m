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
% for
% x_1 >= -3; x_1 <= 3
% x_2 >= -2; x_2 <= 2
%
% Use function fmincon to solve the problem
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and required results.

fun = @(x) (4 - 2.1*x(1)^2/10 + x(1)^4/3)*x(1)^2 + x(1)*x(2) + (- 4 + 4*x(2)^2)*x(2)^2;
x0 = [-3,-2];
x = fmincon(fun,x0,[],[],[],[],[-3 -2],[3 2]);
disp("Problem 1:");
disp(strcat("Find the minimum value starting from the point [",num2str(x0(1)),",",num2str(x0(2)),"]"));
disp(strcat("x_1 = ",num2str(x(1))));
disp(strcat("x_2 = ",num2str(x(2))));
disp(" ");

% objective = (4 - 2.1 * x1^2 + x1^4/3) * x1^2 + x1 * x2 + (-4 + 4 * x2^2) * x2^2;
% x0 = [-2,5 2.5]';
% disp(['Initial Objective: ' num2str(objective(x0))])
% 
% function [c, ceq] = nlcon(x)
%     c = 2;
%     ceq = 4;
% end
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% lb = [-3 -2]';
% ub = [3 2]';
% nonlcon = @nlcon;
% 
% [x, fval, ef, output, lambda] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, nonlcon);
% 
% disp(x)
% disp(['Final Objective: ' num2str(objective(x))])

%% ************************************************************************
% Problem 2:
% Using function fminsearch minimize Branin’s function:
%
%                      2             2
% f(x) = a (x_2 - b x_1  + c x_1 - r)  + s (1 - t) cos(x_1) + s
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

fcnMin = @(x) a*(x(2) - b*x(1)^2 + c*x(1) - r)^2 + s*(1 - t)*cos(x(1)) + s;
x_guess = [15 10];
xmin = fminsearch(fcnMin,x_guess);
disp("Problem 2:");
disp(strcat("Minimize the function with starting point [",num2str(x_guess(1)),",",num2str(x_guess(2)),"]"));
disp(strcat("x_1 = ",num2str(xmin(1))));
disp(strcat("x_2 = ",num2str(xmin(2))));
disp(" ");

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
