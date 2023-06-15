function [] = myplot(alpha, beta,scale)

data = [-1,-1,1,1];

loc = [-beta, -alpha, alpha, beta].*scale;

prange = scale*linspace(-alpha-beta,alpha+beta,20);

% Compute polynomial coefficients
% by solving linear system
% a0 + a1*x + a2*x^2 + a3*x3

C = [ones(4,1),loc',loc'.^2,loc'.^3];

coeff = C\data';

f = @(x) [ones(size(x)),x,x.^2,x.^3] * coeff;

figure
plot(prange,f(prange'));

% Using lagrange polynomial interpolation

l = @(x) data(1)*((x-loc(2)).*(x-loc(3)).*(x-loc(4)))./((loc(1)-loc(2)).*(loc(1)-loc(3)).*(loc(1)-loc(4))) +...
         data(2)*((x-loc(3)).*(x-loc(4)).*(x-loc(1)))./((loc(2)-loc(3)).*(loc(2)-loc(4)).*(loc(2)-loc(1))) +...
         data(3)*((x-loc(4)).*(x-loc(1)).*(x-loc(2)))./((loc(3)-loc(4)).*(loc(3)-loc(1)).*(loc(3)-loc(2))) +...
         data(4)*((x-loc(1)).*(x-loc(2)).*(x-loc(3)))./((loc(4)-loc(1)).*(loc(4)-loc(2)).*(loc(4)-loc(3)));

figure
plot(prange,l(prange))

derivative = coeff(2)

derivative2 = ((data(3)-data(2))*beta^2/2/alpha - (data(4)-data(1))*alpha^2/2/beta)/(beta^2-alpha^2)/scale

