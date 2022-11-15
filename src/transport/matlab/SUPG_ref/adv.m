clc;clear

%N = [4,8,16,64];
N = [64];
M = N+1;
% Change a here
a = 200; k = 1;

fexact = @(x) 24/a*(k/a)^3 + 24/a*(k/a)^2*x + 12/a*(k/a)*x.^2 + 4/a*x.^3 + ...
              (1 - 24/a*(k/a)^3 - 24/a*(k/a)^2 - 12/a*(k/a) - 4/a)/exp(a/k) * exp(a/k.*x);

% Define Gauss quadrature points and weights
gp = [-sqrt(3/5),0.0,sqrt(3/5)];
gw = [5/9,8/9,5/9];

map = @(x,x1,x2) (x1+x2)/2 + x*(x2-x1)/2;

% Define Peclet number
Pe = a/k;

% Set Tau to be 0 for Galerkin formulation
Tau = @(h) h/(2*a)*(coth(Pe) - 1.0/Pe);
%Tau = @(h) 0;

C = 50;

% Setup Nitsche parameters
Sk = 1;
Vk = @(h) C*1/h*(1+abs(Sk));

err = zeros(4,1);

for j = 1:size(N,2)

    h = 1.0/N(j);
    x = [0,(1:N(j))*h];
    % Define SUPG stabilization term
    %tau = h/(2*a)*(coth(Pe) - 1.0/Pe);
    tau = Tau(h);

    % Create linear system
    A = zeros(M(j),M(j));
    b = zeros(M(j),1); 

    for e = 1:N(j)
        % Gauss quadrature
        for g = 1:3 
            p = map(gp(g),(e-1)*h,e*h);

            A(e,e) = A(e,e)+ gw(g)*basis(e,p,M(j),h,2)*( -a*basis(e,p,M(j),h,1) + k*basis(e,p,M(j),h,2) ) + ...
                     gw(g)*tau*a^2*basis(e,p,M(j),h,2)*basis(e,p,M(j),h,2);

            A(e,e+1) = A(e,e+1)+ gw(g)*basis(e,p,M(j),h,2)*( -a*basis(e+1,p,M(j),h,1) + k*basis(e+1,p,M(j),h,2) ) + ...
                       gw(g)*tau*a^2*basis(e,p,M(j),h,2)*basis(e+1,p,M(j),h,2);

            A(e+1,e) = A(e+1,e)+ gw(g)*basis(e+1,p,M(j),h,2)*( -a*basis(e,p,M(j),h,1) + k*basis(e,p,M(j),h,2) ) + ...
                       gw(g)*tau*a^2*basis(e+1,p,M(j),h,2)*basis(e,p,M(j),h,2);

            A(e+1,e+1) = A(e+1,e+1)+ gw(g)*basis(e+1,p,M(j),h,2)*( -a*basis(e+1,p,M(j),h,1) + k*basis(e+1,p,M(j),h,2) ) + ...
                         gw(g)*tau*a^2*basis(e+1,p,M(j),h,2)*basis(e+1,p,M(j),h,2);

            b(e,1) = b(e,1)+ gw(g)*12*p^2*basis(e,p,M(j),h,1) + gw(g)*12*p^2*tau*a*basis(e,p,M(j),h,2);

            b(e+1,1) = b(e+1,1)+ gw(g)*12*p^2*basis(e+1,p,M(j),h,1) + gw(g)*12*p^2*tau*a*basis(e+1,p,M(j),h,2);
        end
    end
% Assign boundary condition to linear system

% Solve with reduced system
%{
 {A00 = A(1:end-1,1:end-1);
 {A0g = A(1:end-1,end);
 {F = b(1:end-1);
 {newans = A00\(F-A0g*1);
 {newans = [newans;1];
 {
 %}
% Solve for H
%{
 {Ag0 = A(end,1:end-1);
 {Fg = b(end);
 {Agg = A(end,end);
 {Mg = 1;
 {H = (Ag0*newans(1:end-1) - Fg + Agg)/Mg;
 {
 %}

A = A*h/2;
b = b*h/2;

% Solve with Nitsche method
endg = M(j);
w1 = basis(endg,1.0,M(j),h,1);
w2 = basis(endg-1,1.0,M(j),h,1);
dw2 = basis(endg-1,1.0,M(j),h,2);
dw1 = basis(endg,1.0,M(j),h,2);

nitsche1 = w1*(a*w1-k*dw1) - Sk*k*dw1*w1 + Vk(h)*w1*w1;
nitsche2 = Sk*k*dw1 - Vk(h)*w1;
nitsche3 = w1*(a*w2-k*dw2) - Sk*k*dw1*w2 + Vk(h)*w1*w2;
nitsche4 = w2*(a*w1-k*dw1) - Sk*k*dw2*w1 + Vk(h)*w2*w1;
nitsche6 = Sk*k*dw2 - Vk(h)*w2;

A(end,end) = A(end,end) + nitsche1;
A(end,end-1) = A(end,end-1) + nitsche3;
A(end-1,end) = A(end-1,end) + nitsche4;

b(end) = b(end) - nitsche2;
b(end-1) = b(end-1) -nitsche6;

sol = A\b;

% Plot
figure
hold on
plot(x,sol,'-o')
plot(linspace(0,1,100),fexact(linspace(0,1,100)),'-')
[t,s] = title(['The Peclet number is ',num2str(Pe), ' ,h^{-1} is ',num2str(N(j)), ' , C is ', num2str(C)]);
s.FontAngle = 'italic';
legend({'FE solution','Exact solution'},'Location','northwest');
hold off

% Calculation of L2 error
for e = 1:N(j)
        % Gauss quadrature
        for g = 1:3 
            p = map(gp(g),(e-1)*h,e*h);

            err(j) = err(j)+gw(g)*(sol(e)*basis(e,p,M(j),h,1)+ sol(e+1)*basis(e+1,p,M(j),h,1) - fexact(p))^2; 

        end
end
err(j) = err(j)*h/2;

%{
 {figure
 {hold on
 {for q=1:N(j)+1
 {        plot(x,basis(q,x,M(j),h,1));
 {end
 {title("Linear basis function");
 %}
end

%{
 {figure
 {loglog(1./N,err)
 {title(['L_2 error ,C= ',num2str(C)]);
 {xlabel('h');
 {ylabel('error');
 %}
