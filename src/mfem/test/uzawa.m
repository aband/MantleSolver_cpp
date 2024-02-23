function sol = uzawa(A,B,C,F,G,tol,tau,maxIter)

r = 1.0;
iter = 0;

x = zeros(size(F));
y = zeros(size(G));

while (r>tol) && (iter < maxIter)

    iter + 1;

    tmp1 = A\(F - (A*x - B'*y));

    x = x + tmp1;

    tmp2 = -B*x - C*y + G;
    %tmp2 = tmp2 - mean(tmp2);

    y = y + tau*tmp2;

    iter = iter +1;

    r = norm(tmp1) + norm(tmp2)

end

sol = [x;y];

iter
