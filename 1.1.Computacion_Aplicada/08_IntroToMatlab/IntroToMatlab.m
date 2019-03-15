A = [2 -3 4; -1 5 2; 3 1 -1];
B = [3 1 2; -1 -1 2; -3 1 1];
C = [0.1 0.2; 1 0.1]
b = [-1; 12; -10];
x = inv(A)*b;

"Answer to problem 1:"
A+B
"Answer to problem 2:"
b'*B
"Answer to problem 3:"
inv(A)
"Answer to problem 4:"
x
"Answer to problem 5:"
dot(A,B)
"Answer to problem 6:"
A*B
"Answer to problem 7:"
exp(C)
expm(C)
"Answer to problem 8:"
log(C)
logm(C)
"Answer to problem 9:"
sqrt(A)
sqrtm(A)