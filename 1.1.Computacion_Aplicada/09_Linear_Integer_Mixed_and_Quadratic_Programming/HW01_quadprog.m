% Problem 2 of Hw 1

warning('off')
clc;

load('Tiingo_data.mat')
D = data;

a1 = D(:,1);

for k = 2:10
    r1(k) = (a1(k)-a1(k-1))/a1(k);
end
%%
r1 = (a1(2:10)-a1(1:9))./a1(1:9);
%%
R = (D(2:end,:)-D(1:end-1,1))./D(1:end-1,:);
mu = mean (R)';
%%
C = cov(R);
%%
N = 9;
k = 1000;
f = -mu;
H = 2*C*k;
Aeq = ones(1,N);
beq = 1;
lb = zeros(1,N);
ub = ones(1,N);
x = quadprog(H,f,[],[],Aeq,beq,lb,ub);
disp(x)

%% Expected Return
er = x'*mu;
disp(er)

%% Expected Variance
ev = x'*C*x;
disp(ev)