%% ************************************************************************
% * AUTHOR(S) :
% *     Bruno González Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW03.m
% *
% * DESCRIPTION :
% *     Computación Aplicada (Ene 19 Gpo 1)
% *     Final Exam
% *
% * NOTES :
% *     In submitting the solution to this final exam, We Bruno González
% *     Soria and Antonio Osamu Katagiri Tanaka affirm our awareness of the
% *     standards of the Tecnológico de Monterrey Ethics Code.
% *
% * START DATE :
% *     02 May 2019
%% ************************************************************************
% This script should start with the command rng(31416), and should not
% contain any other call that initializes the state of the random number
% generator. 

close all, clear all, clc, format compact
rng(31416)

%% ************************************************************************
% Problem 1: OPTIMIZATION
% Consider the following function:
% 
%                  ____6___                     2
%                  \      /             18  / ixi \
%         f(X) =    \        sin(xi) sin   | ----- |
%                   /                       \  Pi /
%                  /_______\
%                     i=1
% 
% where 0 < xi < 5.
% 
% Maximize function f () using the Nelder-Mead algorithm (fminsearch) and
% simulated annealing (simulannealbnd). Modify whatever parameters you deem
% necessary to produce a good performance of these algorithms, regardless
% of the state of the random number generator. Use randomly generated
% initial point in the valid range of x.
%% ************************************************************************
%   a) Implement f (x) as a MATLAB function.
i = 1:6;
f = @(x) -fx(x);

%% ************************************************************************
%   b) Give your best solution found (optimal x and evaluation of x) for
%   each algorithm.

% NelderMeade
x0 = rand([1 6])*5;
options = optimset('Display', 'off', 'MaxFunEvals', 10000);
disp("The optimal value of x usning NelderMeade (fminsearch) method is:")
[x,fval,exitflag,output] = fminsearch(f,x0,options)

% Simulated Annealing
disp("The optimal value of x usning Simulated Annealing (simulannealbnd) method is:")
lb = zeros([1 6]);
ub = ones([1 6])*5;
[x,fval,exitflag,output] = simulannealbnd(f,x0,lb,ub,options)

%% ************************************************************************
%   c) Which of these two algorithms has a better expected performance on
%   this problem when varying the initial point(s)? Justify your answer.

% From this two algorithms, fminsearch has a better expected performance on
% this problem. The reason is that the fval (objective function value at
% the solution) obtained is larger than the one obtained in simulannealbnd,
% thus closer to a maximum in the function. Additionally, unlike other
% solvers, fminsearch stops when it satisfies both TolFun and TolX.

%% ************************************************************************
%   DEFINED FUNCTIONS:
%   Problem 1 a)
function fcn = fx (x)
suma=0;
    for i = 1:6
        newterm = sin(x(i))*sin((i*x(i))^2/pi)^18;
        suma = suma + newterm;
    end
    fcn = suma;
end

