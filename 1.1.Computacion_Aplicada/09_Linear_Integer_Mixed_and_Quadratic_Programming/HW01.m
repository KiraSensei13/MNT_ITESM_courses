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
% *     Homework on Linear, mixed, and quadratic programming 
% *
% * NOTES :
% *     
% *
% * START DATE :
% *     04 Apr 2019
% ************************************************************************

clc;
clear all;

%% ************************************************************************
% Problem 1:
% Solve the Following transportation problem:
s = [37.6; 40.4; 44.5];
d = [20 30 30 40];
C = [41 27 28 24; 40 29 100 23; 37 30 27 21];

% transportation costs tableau
%             | destination     | supply
%             |  1   2   3   4  | 
% ------------+-----------------+--------
%    source 1 |  41  27  28  24 | 37.6
%           2 |  40  29 100  23 | 40.4
%           3 |  37  30  27  21 | 44.5
% ------------+-----------------+--------
%      demand |  20  30  30  40 |

% Assume that only integer units can be transported.
% Assume that only multiples of 2 units can be transported.
% Upload to Blackboard a pdf file with a MATLAB script that solves both
% cases.
% The pdf should also include the solution matrix A for each case.

% Only multiples of 2 units can be transported, so let's assume each unit
% contains 2 objects. Let's simulate this by divition by 2.
s = s/2;
d = d/2;
C = C/2;

% Only integer units can be transported, so round toward negative infinity.
s = floor(s);

% The following is based on:
% http://web.tecnico.ulisboa.pt/mcasquilho/compute/_linpro/TaylorB_module_b.pdf
f = [C(1,:) C(2,:) C(3,:)]';
n = length(s);
m = length(d);
f = reshape(C',n*m,1);
A = zeros(n,n*m);
%
for i=1:n
    A(i,1+(i-1)*4:i*4) = 1;
end
%
b = s;
Aeq = zeros(n,n*m);
%
for j=1:m
    Aeq(j,j:m:n*m) = 1;
end
%
beq = d;
LB = zeros(n*m,1);
UB = Inf(n*m,1);
%
x = linprog(f,A,b,Aeq,beq,LB,UB);
%
assig = reshape(x,m,n)';
total_cost = sum(sum(assig.*C));

% Up to this point calculations have been made by units of 2 objects, let's
% multiply the result by 2 to get the transportation assignments and const
% by product
assig = assig*2;
total_cost = total_cost*2;

% Print the calculations
disp(assig);
disp(total_cost);

clear all;

%% ************************************************************************
% Problem 2:
% Using function Investopedia, download the adjusted closing prices for the
% following DJI stocks from January 1, 2018 through January 1, 2019.
%
% KO Coca-Cola            DIS Disney
% PG Procter & Gamble     MCD McDonald’s
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

% Expected returns / mean of the returns per security
muR  = mean(Returns);

% Returns variance per security / risk
C = cov(Returns); % std(R)^2 variances are in the main diagonal ...
sigmaR = diag(C); % extract the main diagonal of C

% Plot their expected return versus variance.
figure(1)
scatter(sigmaR,muR,1)
title("Pareto Front")
xlabel("risk") 
ylabel("return") 
fields = fieldnames(res);
nAssets = length(fields); % number of assets
for i = 1:nAssets
	text(sigmaR(i),muR(i),strcat( "  ", fields(i)) );
end

% Quadratic Programming Portfolio Optimization

