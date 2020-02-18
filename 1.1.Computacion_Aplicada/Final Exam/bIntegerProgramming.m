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
% Problem : INTEGER PROGRAMMING
% 
% An airline company is considering the purchase of new long-, medium-, and
% short-range jet passenger airplanes. The purchase price is $33.5M for
% each long-range plane, $25M for each medium-range plane, and $17.5M for
% each short-range plane. The board of directors has authorized a maximum
% of $750M for these purchases. Regardless of which planes are purchased,
% air travel of all distances is expected to be sufficiently large enough
% so that these planes would be utilized at essentially maximum capacity.
% It is estimated that the net annual profit (after subtracting capital
% recovery costs) would be $2.1M per long-range plane, $1.5M per
% medium-range plane, and $1.15M per short-range plane.
% 
% Enough trained pilots are available to the company to crew 30 new
% airplanes. If only short-range planes were purchased, facilities would be
% able to handle 40 new planes. However, each medium-plane is equivalent to
% 1+1/3 short-range planes, and each long-range plane is equivalent to
% 1+2/3 short-range planes in terms of their use of maintenance facilities.
% Using the preceding data, management wishes to know how many planes of
% each type should be purchased to maximize profit.

%% ************************************************************************
% a) Formulate the problem as an integer programming problem.
%
%    Let L be the number of long-range jets to buy
%    Let M be the number of medium-range jets to buy
%    Let S be the number of short-range jets to buy
%    And let P the profit
%    
%    Maximise P = 2.1*L + 1.5*M + 1.15*S
%
%    Subject to: 33.5*L + 25*M + 17.55*S <= 750
%                L + M + S <= 30
%                1.67*L + 1.33*M + S <= 40
%                L >= 0, M >= 0, S >= 0
%                L, M, S are integers
%
% b) Use intlinprog to find the solution (number of planes of each type and
%    maximum profit)

l = optimvar('L','Type','integer');
m = optimvar('M','Type','integer');
s = optimvar('S','Type','integer');
prob = optimproblem('ObjectiveSense','maximize');
prob.Objective = 2.1*l + 1.5*m + 1.15*s;
prob.Constraints.cons1 = 33.5*l + 25*m + 17.55*s <= 750;
prob.Constraints.cons2 = l + m + s <= 30;
prob.Constraints.cons3 = (5/3)*l + (4/3)*m + s <= 40;
prob.Constraints.cons4 = l >= 0;
prob.Constraints.cons5 = m >= 0;
prob.Constraints.cons6 = s >= 0;

options = optimoptions('intlinprog','Display','off');

sol = solve(prob);
disp(sol);
