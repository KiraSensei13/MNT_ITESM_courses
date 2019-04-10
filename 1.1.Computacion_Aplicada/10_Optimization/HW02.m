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


