% ************************************************************************
% * AUTHOR(S) :
% *     Bruno González Soria          (A01169284)
% *     Antonio Osamu Katagiri Tanaka (A01212611)
% *
% * FILENAME :
% *     HW02.m
% *
% * DESCRIPTION :
% *     Computaciónn Aplicada (Ene 19 Gpo 1)
% *     Homework on Optimization 
% *
% * NOTES :
% *     
% *
% * START DATE :
% *     11 Apr 2019
% ************************************************************************

clc;
% clear all;

% syms -> Short-cut for constructing symbolic variables.
% diff -> Differentiates a symbolic expression.
% solve -> Symbolic solution of algebraic equations.
% subs -> Symbolic substitution
% fminunc -> Unconstrained minimization. Finds a local minimum of a function of several variables.
% fminbnd -> Single-variable bounded nonlinear function minimization. Finds a local minimal in the
% interval -> xmin < x < xmax.
% fmincon -> Finds a constrained minimum of a function of several variables.
% fminsearch -> Multidimensional unconstrained nonlinear minimization (Nelder-Mead).

%% ************************************************************************
% Problem 1:
% Solve the following problem using the optimization toolbox:
% min f(x) = (4 - 2.1 x^2_1 + x^4_1/3) x^2_1 + x_1 x_2 + (-4 + 4x^2_2) x^2_2
% for x_1 [-3; 3] and x_2 [-2; 2].
% 1 Use function fmincon to solve the problem.
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and required results.

objective = (4 - 2.1 * x1^2 + x1^4/3) * x1^2 + x1 * x2 + (-4 + 4 * x2^2) * x2^2;
x0 = [rand*15-5 rand*15]';
xmin = fmincon(objective,x0,[],[],[],[],[-3 3]',[-2 2]');
print xmin

%% ************************************************************************
% Problem 2:
% Using function fminsearch minimize Branin’s function:
% f(x) = a(x^2 - bx^2_1 + cx_1 - r)^2 + s(1 - t) cos(x_1) + s
% where
% a = 1             b = 5.1 / 4Pi^2
% c = 5 / Pi        r = 6
% s = 10            t = 1 / 8Pi
% for x_1 [-5; 10] y x_2 [0; 15].
% Upload to Blackboard a pdf file that contains a MATLAB script, any MATLAB
% functions that you implemented, and required results.

function y = branin(xx)
% Branin's function
x1 = xx(1);
x2 = xx(2);
t = 1 / (8*pi);
s = 10;
r = 6;
c = 5/pi;
b = 5.1 / (4*pi^2);
a = 1;
term1 = a * (x2 - b*x1.^2 + c*x1 - r).^2;
term2 = s*(1-t)*cos(x1);
y = term1 + term2 + s;
end

function y = branin2(x1,x2)
% Branin's function (used for plotting)
t = 1 / (8*pi);
s = 10;
r = 6;
c = 5/pi;
b = 5.1 / (4*pi^2);
a = 1;
term1 = a * (x2 - b*x1.^2 + c*x1 - r).^2;
term2 = s*(1-t)*cos(x1);
y = term1 + term2 + s;
end

%% Solution using fminsearch
x0 = [rand*15-5 rand*15]';
xmin = fmincon(@branin,x0,[],[],[],[],[-5 0]',[10 15]');
branin(xmin)

%% Branin's function graph
[X1,X2] = meshgrid(linspace(-5,10),linspace(0,15));
F = branin2(X1,X2);
surf(X1,X2,F)
shading interp
xlabel('x_1')
ylabel('x_2')