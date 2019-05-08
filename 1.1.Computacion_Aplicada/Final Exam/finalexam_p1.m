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
% *     In submitting the solution to this final exam, I (we) ?your name(s)?
% *     affirm my (our) awareness of the standards of the Tecnológico de
% *     Monterrey Ethics Code.
% *
% * START DATE :
% *     02 May 2019
%% ************************************************************************
% This script should start with the command rng(31416), and should not
% contain any other call that initializes the state of the random number
% generator. 

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



%% ************************************************************************
%   b) Give your best solution found (optimal x and evaluation of x) for
%   each algorithm.

%% ************************************************************************
%   c) Which of these two algorithms has a better expected performance on
%   this problem when varying the initial point(s)? Justify your answer.

