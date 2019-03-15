clc

A = [ 1 -3  4;
     -1  5  2;
      3  1 -1];

B = [ 3  1 2;
     -1 -1 2;
     -3  1 1];
 
b = [ -1;
      12;
     -10];
 
C = [0.1 0.2;
     1   0.1];

A + B
b.' * B
inv(A)
% Solve A*x = b for x
x = A ./ b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dot(A,B)
cross(A,B)

exp(C)
expm(C)

log(C)
logm(C)

sqrt(C)
sqrtm(C)

