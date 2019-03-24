% clc
% 
% A = [ 1 -3  4;
%      -1  5  2;
%       3  1 -1];
% 
% B = [ 3  1 2;
%      -1 -1 2;
%      -3  1 1];
%  
% b = [ -1;
%       12;
%      -10];
%  
% C = [0.1 0.2;
%      1   0.1];
% 
% A + B
% b.' * B
% inv(A)
% % Solve A*x = b for x
% x = A ./ b
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% dot(A,B)
% cross(A,B)
% 
% exp(C)
% expm(C)
% 
% log(C)
% logm(C)
% 
% sqrt(C)
% sqrtm(C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% An anominous function ...
f = @(x) exp(-x./4).*sin(x + 0.2) + 0.1.*x;
f(1)
f([1 2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linear programming example
f = @(x_1,x_2) 3 - 2*x_1 + 1.2*x_2;

f = [-2;1.2];
A = [2 -1];
b = 2.5;
LB = [0; 0];
linprog(f,A,b,[],[],LB)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Transportation Problem

s = [75; 75; 45];






