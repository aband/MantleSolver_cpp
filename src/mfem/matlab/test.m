clc; clear

x = linspace(0,1,10);
y = linspace(0,1,10);

[X,Y] = meshgrid(x,y);

% Interpolant f(x,y) = a0 + a10 x + a01 y + a11 xy ...
% Sample points (0,0), (0.5,0.4) (0.3,0.6)
% Interpo values 1, 1, 1

shift = @(a,b,x) 0.5*(a+b) + 0.5*(b-a)*x
