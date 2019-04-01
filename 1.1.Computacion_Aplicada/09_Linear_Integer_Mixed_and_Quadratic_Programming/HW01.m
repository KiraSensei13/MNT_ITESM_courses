% ************************************************************************
% * AUTHOR(S) :
% *     Bruno Gonz�lez Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW01.m
% *
% * DESCRIPTION :
% *     Computaci�nn Aplicada (Ene 19 Gpo 1)
% *     Homework on Linear, mixed, and quadratic programming 
% *
% * NOTES :
% *     
% *
% * START DATE :
% *     04 Apr 2019
% ************************************************************************
% Problem 1:
% Solve the Following transportation problem:
s = [37.6; 40.4; 44.5];
d = [20 30; 30; 40];
C = [41 27 28 24; 40 29 100 23; 37 30 27 21];
% Assume that only integer units can be transported.
% Assume that only multiples of 2 units can be transported.
% Upload to Blackboard a pdf file with a MATLAB script that solves both cases.
% The pdf should also include the solution matrix A for each case.


% ************************************************************************
% Problem 2:
% Using function Investopedia, download the adjusted closing prices for the
% following DJI stocks from January 1, 2018 through January 1, 2019.
%
% KO Coca-Cola            DIS Disney
% PG Procter & Gamble     MCD McDonald�s
% PFE Pfizer              WMT WalMart
% MRK Merck               V Visa
% VZ Verizon
%
% Obtain the daily returns for these securities. Plot their expected return
% versus variance.
% Using command quadprog, obtain the optimal portfolios for
% k = 1; 2; 2:5; 4; 7; 9; 11; 20; 50; 1000
% For each value of k, report optimal weights, expected return and variance.
% Plot this data on the previous graph.
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and the required plots and results.