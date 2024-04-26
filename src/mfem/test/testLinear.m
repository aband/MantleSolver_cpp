% Test 

%function sol = testLinear(tol, tau, maxIter)
clc; clear

readMat;
format long;

% Form saddle point linear system passed into uzawa solver

tol     = 10e-9;
tau     = 1;
maxIter = 500;

lA = As;

M = size(Bs,2);
N = size(Bs,1);

lB = [Bs;zeros(N,M)];

% Create schur complementa
S = Bd*inv(Ad)*Bd' + Cd;

lC = [Cs,K;K,S];

sol = uzawa(lA, lB, lC, gs1', G, tol, tau, maxIter);
