function p3 = p3reconst(N)
format long
% Get mesh
x = linspace(0,1,4);
y = linspace(0,1,4);

%{
 {x = [0,0.5,1.0,1.5;0,0.56601,1.0,1.5;0,0.5,1.0,1.5];
 {y = [0,0,0,0; 0.5,0.56601,0.5,0.5;1.0,1.0,1.0,1.0];
 {
 %}
% Normalize

% Center reference point (0.5,0.5)

a0 = @(x,y) 0*x+1;
a1 = @(x,y) x;
a2 = @(x,y) x.^2;
a3 = @(x,y) y;
a4 = @(x,y) y.*x;
a5 = @(x,y) y.*(x.^2);
a6 = @(x,y) y.^2;
a7 = @(x,y) (y.^2).*x;
a8 = @(x,y) (y.^2).*(x.^2);



%fun = @(x,y) sin(x).*sin(y);
fun = @(x,y) x.^2;

exact   = zeros(3,3);
cellave = zeros(3,3);

x = x/N;
y = y/N;

% cell coordinates
cx = (x(1:end-1) + x(2:end))/2;
cy = (y(1:end-1) + y(2:end))/2;


area = integral2(a0,x(1),x(2),y(1),y(2));

for j=1:3
    for i = 1:3
        exact(j,i)   = fun(cx(i),cy(j));
        cellave(j,i) = integral2(fun,x(i),x(i+1),y(j),y(j+1))/area;
    end
end

x = (x - 0.5)*3;
y = (y - 0.5)*3;

% Calculate region integrals

area = integral2(a0,x(1),x(2),y(1),y(2));

% exact solution and cell averaged values
% Calculate region integrals

M = zeros(9,9);


sol = zeros(9,9);

for k = 1:9
    for j=1:3
        for i=1:3
            % Rescale coordinates

            ksi0 = x(i);
            ksi1 = x(i+1);
            eta0 = y(j);
            eta1 = y(j+1);

            M((j-1)*3+i,1)  = integral2(a0,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,2)  = integral2(a1,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,3)  = integral2(a2,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,4)  = integral2(a3,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,5)  = integral2(a4,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,6)  = integral2(a5,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,7)  = integral2(a6,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,8)  = integral2(a7,ksi0,ksi1,eta0,eta1);
            M((j-1)*3+i,9)  = integral2(a8,ksi0,ksi1,eta0,eta1);
		  end
	 end
M = M/area

B = zeros(9,1);
B(k) = 1;

sol(k,:) = M\B;
end

sol
sol(1,:)

exact

cellave

% Reconstruction and sampling
xr = x(1)
yr = y(1)
rec = 0;

for j=1:3
		  for i=1:3
rec = rec + cellave(j,i)*(sol((j-1)*3+i,1)*a0(xr,yr) + ...
                          sol((j-1)*3+i,2)*a1(xr,yr) + ...
                          sol((j-1)*3+i,3)*a2(xr,yr) + ...
                          sol((j-1)*3+i,4)*a3(xr,yr) + ...
                          sol((j-1)*3+i,5)*a4(xr,yr) + ...
                          sol((j-1)*3+i,6)*a5(xr,yr) + ...
                          sol((j-1)*3+i,7)*a6(xr,yr) + ...
                          sol((j-1)*3+i,8)*a7(xr,yr) + ...
                          sol((j-1)*3+i,9)*a8(xr,yr) );
		  end
end

res =  rec - fun(xr,yr)

end
