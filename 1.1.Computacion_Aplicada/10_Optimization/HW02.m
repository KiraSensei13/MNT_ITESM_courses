% ************************************************************************
% * AUTHOR(S) :
% *     Bruno Gonzï¿½lez Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW01.m
% *
% * DESCRIPTION :
% *     Computación Aplicada (Ene 19 Gpo 1)
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
%            /                   4 \
%            |            2   x_1  |    2                         2     2
% min f(x) = | 4 - 2.1 x_1 +  ---- | x_1  + x_1 x_2 + (- 4 + 4 x_2 ) x_2
%  x         \                  3  /
%
% for
% x_1 >= -3; x_1 <= 3
% x_2 >= -2; x_2 <= 2
%
% Use function fmincon to solve the problem
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and required results.

fun = @(x) (4 - 2.1*x(1)^2 + x(1)^4/3)*x(1)^2 + x(1)*x(2) + (- 4 + 4*x(2)^2)*x(2)^2;
x0 = [-3,-2];
x = fmincon(fun,x0,[],[],[],[],[-3 -2],[3 2]);
disp("Problem 1:");
disp(strcat("Find the minimum value starting from the point [",num2str(x0(1)),",",num2str(x0(2)),"]"));
disp(strcat("x_1 = ",num2str(x(1))));
disp(strcat("x_2 = ",num2str(x(2))));
disp(" ");

%% ************************************************************************
% Problem 2:
% Using function fminsearch minimize Braninï¿½s function:
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
x_guess = [rand*15-5 rand*15]';
xmin = fminsearch(fcnMin,x_guess);
disp("Problem 2:");
disp(strcat("Minimize the function with starting point [",num2str(x_guess(1)),",",num2str(x_guess(2)),"]"));
disp(strcat("x_1 = ",num2str(xmin(1))));
disp(strcat("x_2 = ",num2str(xmin(2))));
disp(" ");

