% 1D code utilizing WENO reconstruction for FV scheme
% A boundary value problem
% Finite volume method cannot solve pure Neumann problem prescribing total flux

clc; clear

% Define cell number
N = 64;

% Define vertex number
M = N + 1;

% Define right hand side

x = linspace(0,1,M);

f = 4*(x(2:end).^3 - x(1:end-1).^3);

% Define convection and diffusion coefficients
% A convection dominated convection-diffusion problem
a = 200; k = 1;

% A WENO (3,2) reconstruction
% Define reconstruction stencls

stencil32 = [[-1,1]; [-1,0]; [0,1];
linWgt32  = [3,1,1];


