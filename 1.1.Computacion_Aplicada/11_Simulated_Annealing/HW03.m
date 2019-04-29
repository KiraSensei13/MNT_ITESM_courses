% ************************************************************************
% * AUTHOR(S) :
% *     Bruno Gonzï¿½lez Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW03.m
% *
% * DESCRIPTION :
% *     Computación Aplicada (Ene 19 Gpo 1)
% *     Homework on Simulated Annealing
% *
% * NOTES :
% *     
% *
% * START DATE :
% *     25 Apr 2019
% ************************************************************************

warning('off')
clc;
clear all;
close all;

%% ************************************************************************
% Problem 1:
% Using simulated annealing in the global optimization toolbox, solve the
% traveling salesman problem for 70 cities with coordinates generated in the
% following way:
N = 70;
rng(123);
coordinates = rand(N,2);
% Upload to Blackboard a pdf file that contains the folowing:
%   A script of your solution
%   A plot of the best route found

distances = TSPtable(coordinates);
objFcnTSP(distances);
route0 = randperm(N)';
TSPplot(route0,coordinates,'b',2)
xlabel('x')
ylabel('y')
title(sprintf('TSP?with?%d?cities',N))
options = optimoptions(@simulannealbnd,'DataType','custom',...'AnnealingFcn',@TSPinversion);
r0 = randperm(N);
TSPplot(r0,coordinates,'-')
title(sprintf('initial random route, cost=%f',objFcnTSP(r0)))
r = simulannealbnd(@objFcnTSP,r0,[],[],options);
TSPplot(r,coordinates,'-')
title(sprintf('cost=%f',objFcnTSP(r)))

%% ************************************************************************
% Problem 2:
% Using function Tiiago, download the adjusted closing prices for the following DJI stocks from
% January 1, 2018 through January 1, 2019.

% KO Coca-Cola                     DIS Disney
% PG Procter & Gamble              MCD McDonald’s
% PFE Pfizer                       WMT WalMart
% MRK Merck                        V Visa
% VZ Verizon

% Weights must be in the range [0:1; 0:5].
% There must be at least 3 and no more than 5 securities in the portfolio.
% Obtain the daily returns for these securities. Plot their expected return versus variance.
% Using command simulated annealing, obtain the optimal portfolios for
% k = 1; 2; 2:5; 4; 7; 9; 11; 20; 50; 1000
% For each value of k, report optimal weights, expected return and variance.
% Plot this data on the previous graph.
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB functions that you
% implemented, and the required plots and results.

load('Tiingo_data.mat')
% Let's convert the closing prices into returns

% CLOSING PRICES
CP_KO  = res.KO;
CP_PG  = res.PG;
CP_PFE = res.PFE;
CP_MRK = res.MRK;
CP_VZ  = res.VZ;
CP_DIS = res.DIS;
CP_MCD = res.MCD;
CP_WMT = res.WMT;
CP_V   = res.V;

% RETURNS
R_KO  = (CP_KO(2:end)-CP_KO(1:end-1))./CP_KO(1:end-1);
R_PG  = (CP_PG(2:end)-CP_PG(1:end-1))./CP_PG(1:end-1);
R_PFE = (CP_PFE(2:end)-CP_PFE(1:end-1))./CP_PFE(1:end-1);
R_MRK = (CP_MRK(2:end)-CP_MRK(1:end-1))./CP_MRK(1:end-1);
R_VZ  = (CP_VZ(2:end)-CP_VZ(1:end-1))./CP_VZ(1:end-1);
R_DIS = (CP_DIS(2:end)-CP_DIS(1:end-1))./CP_DIS(1:end-1);
R_MCD = (CP_MCD(2:end)-CP_MCD(1:end-1))./CP_MCD(1:end-1);
R_WMT = (CP_WMT(2:end)-CP_WMT(1:end-1))./CP_WMT(1:end-1);
R_V   = (CP_V(2:end)-CP_V(1:end-1))./CP_V(1:end-1);

Returns = [R_KO R_PG R_PFE R_MRK R_VZ R_DIS R_MCD R_WMT R_V];