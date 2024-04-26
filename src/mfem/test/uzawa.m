function sol = uzawa(A,B,C,F,G,tol,tau,maxIter)

r = 1.0;
iter = 0;

x = zeros(size(F));
y = zeros(size(G));

S = B*inv(A)*B';

%S = B*B';

while (r>tol) && (iter < maxIter)

    iter + 1;

    tmp1 = A\(F - (A*x - B'*y));

    x = x + tmp1;

    tmp2 = -B*x - C*y + G;
    %tmp2 = tmp2 - mean(tmp2);

    % =====================================
    % MINRES
    [tmp2, minr] = minres(S,tmp2,zeros(size(tmp2)),5000,1e-10);

    % Unique solution for 
    %tmp2(1) = tmp2(1)/(S(1,1)-S(1,2));
    %tmp2(2) = -tmp2(1);
    % =====================================

    fprintf(' %d : %.10f, %.10f \n',iter, norm(tmp1), norm(tmp2));

    y = y + tau*tmp2;

    iter = iter +1;

    r = norm(tmp1) + norm(tmp2);

end

sol = [x;y];

iter
