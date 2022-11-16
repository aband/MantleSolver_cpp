% 1D code utilizing WENO reconstruction for FV scheme
% A boundary value problem
% Finite volume method cannot solve pure Neumann problem prescribing total flux
% Solving a sample steady advection-diffusion problem.
% A pure dirichlet problem prescribing dirichlet boundary conditions at both sides

clc; clear

% Define cell number
N = 64;

% Define vertex number
M = N + 1;

% Define right hand side

x = linspace(0,1,M);

h = 1/N;

cell = linspace(0,1,N);

sol = zeros(size(cell));

%f = 4*(x(2:end).^3 - x(1:end-1).^3);

% Define convection and diffusion coefficients
% A convection dominated convection-diffusion problem
a = 200; k = 1;

% A WENO (3,2) reconstruction
% Define reconstruction stencls

stencil32 = [[-1,1]; [-1,0]; [0,1]];
linWgt32  = [3,1,1];

% Degenerate multilevel WENO (3,2) - (2,1) reconstruction at boundary
stencil32L = [[0,1];[0,0]];
linWgt32L  = [2,1];

stencil32R = [[-1,0];[0,0]];
linWgt32R  = [2,1];

% Define exact solution
Pe = a/k;
fexact = @(x) 1/(exp(Pe)-1) * (exp(Pe.*x) - 1);





% Plot
figure
hold on
plot(cell,sol,'-o')
plot(linspace(0,1,100),fexact(linspace(0,1,100)),'-')
[t,s] = title(['The Peclet number is ',num2str(Pe), ' ,h^{-1} is ',num2str(N)]);
s.FontAngle = 'italic';
legend({'FV solution','Exact solution'},'Location','northwest');
hold off
