% Test 

function sol = testLinear(tol, tau, maxIter)

readMat;

% Form saddle point linear system passed into uzawa solver

lA = As;

M = size(Bs,2);
N = size(Bs,1);

lB = [Bs,zeros(N,N)];

% Create schur complementa
S = Bd*inv(Ad)*Bd' - Cd;

lC = [Cs,K;K,S];

sol = uzawa(lA, lB, lC, gs1, [gs2';gd2'], tol, tau, maxIter);
