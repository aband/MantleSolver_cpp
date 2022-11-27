% Test reconstruction of diffustion flux.
% Reconstruction of derivatives.

clc;clear

exact = @(x) sin(x);
exactDerv = @(x) cos(x);

stencil43 = [[-2,1];[-2,0];[-1,1]];
linWgt43  = [4,1,1]; 

stencil32 = [[-1,1];[-1,0];[0,1]];
linWgt32  = [3,1,1]; 

% physical parameters
N = 64;

% Define vertex number
M = N + 1;

% Define right hand side
startP = -1;
endP = 1;
x = linspace(startP,endP,M);
h = (endP-startP)/N;
cell = linspace(startP,endP,N);
sol = zeros(size(cell));

% Define convection and diffusion coefficients
% A convection dominated convection-diffusion problem
a = 200; k = 1;
Pe = a/k;

% Gauss quadrature points and weights
gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
gaussWt = [5/18 8/18 5/18]; %half value!

% Create uBar
uBar = 0;

vertxL = x(1:end-1); 
vertxR = x(2:end);
for g = 1:3
    gPt = vertxL + (vertxR-vertxL)*(gaussPt(g)+1)/2;
    uBar = uBar + gaussWt(g)*sin(gPt);
end

% Biased WENO stencil method for Dirichlet boundary condition
eta_bias = [0,0,0,0];

alpha = 0.5; 
beta  = 1.5; 

% For (3,2) reconstruction, a centerShift -0.5 should be implemented
% For (4,3) reconstruction, a centerShift 0.0 should be implemented otherwise
uExact = sin((x(2)+x(3))/2);
% The first reconstruction on target vertx number 3, vertex centered type reconstruction
% with -0.5 relative coordinate, the actual reconstruction point is 1.5
ru = multiLWENO1D(x,h,uBar,stencil43,linWgt43,3,-0.5,0.0,2);
% The seconf reconstruction on target cell number 2, cell centered type reconstruction
% with 0.0 relative coordinate, the actual reconstruction point is also 1.5
ru = multiLWENO1D(x,h,uBar,stencil32,linWgt32,2,0.0,-0.5,1);
% All three values matched. Test passed.
hatX = [-1.5,-0.5,0.5,1.5];

ru = multiLWENO1D(x,h,uBar,stencil43,linWgt43,3,hatX,0.0,2)

dflux = diffFlux(alpha*h,beta*h,ru)

exact = cos(x(3)) 

Error = exact - dflux

% =========================================================================================
function [flux] = diffFlux(alpha,beta,ru)

    flux = ((ru(3)-ru(2))*beta^2/(2*alpha) - ...
            (ru(4)-ru(1))*alpha^2/(2*beta))/ ...
           (beta^2-alpha^2);

end
