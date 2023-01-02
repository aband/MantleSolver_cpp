function p4res = p4reconst(N)

% Get mesh
x = [0.0,1.0,2.0,3.0,4.0];
y = [0.0,1.0,2.0,3.0,4.0,5.0];

% Normalize
x = (x-x(3))/N;
y = (y-y(3))/N;

x
y

% Center reference point (0.0,0.0)
a0 = @(x,y) 0*x+1;
a1 = @(x,y) x;
a2 = @(x,y) x.^2;
a3 = @(x,y) x.^3;
a4 = @(x,y) y;
a5 = @(x,y) y.*x;
a6 = @(x,y) y.*(x.^2);
a7 = @(x,y) y.*(x.^3);
a8 = @(x,y) y.^2;
a9 = @(x,y) (y.^2).*x;
a10 = @(x,y) (y.^2).*x.^2;
a11 = @(x,y) (y.^2).*x.^3;
a12 = @(x,y) y.^3;
a13 = @(x,y) (y.^3).*x;
a14 = @(x,y) (y.^3).*x.^2;
a15 = @(x,y) (y.^3).*x.^3;
a16 = @(x,y) y.^4;
a17 = @(x,y) (y.^4).*x;
a18 = @(x,y) (y.^4).*x.^2;
a19 = @(x,y) (y.^4).*x.^3;

% Calculate region integrals

area = integral2(a0,x(1),x(2),y(1),y(2));

% cell coordinates
cx = (x(1:end-1) + x(2:end))/2;
cy = (y(1:end-1) + y(2:end))/2;

fun = @(x,y) sin(x).*sin(y);

exact   = zeros(5,4);
cellave = zeros(5,4);

% exact solution and cell averaged values
for j=1:5
    for i = 1:4
        exact(j,i)   = fun(cx(i),cy(j));
        cellave(j,i) = integral2(fun,x(i),x(i+1),y(j),y(j+1))/area;
    end
end

% Calculate coefficients

M = zeros(20,20);

sol = zeros(20,20);

for k = 1:20
for j=1:5
    for i=1:4

        ksi0 = x(i);
        ksi1 = x(i+1);
        eta0 = y(j);
        eta1 = y(j+1);

        M((j-1)*4+i,1)  = integral2(a0,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,2)  = integral2(a1,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,3)  = integral2(a2,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,4)  = integral2(a3,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,5)  = integral2(a4,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,6)  = integral2(a5,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,7)  = integral2(a6,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,8)  = integral2(a7,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,9)  = integral2(a8,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,10) = integral2(a9,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,11) = integral2(a10,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,12) = integral2(a11,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,13) = integral2(a12,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,14) = integral2(a13,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,15) = integral2(a14,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,16) = integral2(a15,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,17) = integral2(a16,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,18) = integral2(a17,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,19) = integral2(a18,ksi0,ksi1,eta0,eta1);
        M((j-1)*4+i,20) = integral2(a19,ksi0,ksi1,eta0,eta1);
	  end
end

M = M/area;

b = zeros(20,1);
b(k) = 1;

sol(k,:) = M\b;

end

format short

M*area

sol

% Reconstruction and sampling
xr = x(1);
yr = y(2);
rec = 0;

for j=1:5
    for i=1:4
        rec = rec + cellave(j,i)*(sol((j-1)*4+i,1)*a0(xr,yr) + ...
                                  sol((j-1)*4+i,2)*a1(xr,yr) + ...
											 sol((j-1)*4+i,3)*a2(xr,yr) + ...
                                  sol((j-1)*4+i,4)*a3(xr,yr) + ...
                                  sol((j-1)*4+i,5)*a4(xr,yr) + ...
                                  sol((j-1)*4+i,6)*a5(xr,yr) + ...
                                  sol((j-1)*4+i,7)*a6(xr,yr) + ...
                                  sol((j-1)*4+i,8)*a7(xr,yr) + ...
                                  sol((j-1)*4+i,9)*a8(xr,yr) + ...
                                  sol((j-1)*4+i,10)*a9(xr,yr) + ...
                                  sol((j-1)*4+i,11)*a10(xr,yr) + ...
                                 sol((j-1)*4+i,12)*a11(xr,yr) + ...
                                  sol((j-1)*4+i,13)*a12(xr,yr) + ...
                                  sol((j-1)*4+i,14)*a13(xr,yr) + ...
                                  sol((j-1)*4+i,15)*a14(xr,yr) + ...
                                  sol((j-1)*4+i,16)*a15(xr,yr) + ...
                                 sol((j-1)*4+i,17)*a16(xr,yr) + ...
                                  sol((j-1)*4+i,18)*a17(xr,yr) + ...
                                  sol((j-1)*4+i,19)*a18(xr,yr) + ...
                                  sol((j-1)*4+i,20)*a19(xr,yr)  );
    end
end

rec;

fun(xr,yr);

p4res = rec - fun(xr,yr);

end
