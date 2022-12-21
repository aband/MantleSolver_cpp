function [] = test(level) 
% test routines creating base polynomial coefficients

stencil = [[-1,1];[0,1];[-1,0]];

grid = ([0.1,0.2,0.4,0.7] - 0.3)/level + 0.3;

%grid = ([0.0,0.2,0.4,0.6] - 0.3)/level + 0.3;

basepolyncoeff = basePolyn(stencil,-0.5,grid);

% Gauss quadrature points and weights
gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
gaussWt = [5/18 8/18 5/18]; %half value!

% Create uBar
uBar = 0;

fun = @(x) sin(x);

vertxL = grid(1:end-1); 
vertxR = grid(2:end);
for g = 1:3
    gPt = vertxL + (vertxR-vertxL)*(gaussPt(g)+1)/2;
    uBar = uBar + gaussWt(g)*fun(gPt);
end

integral(fun,grid(1),grid(2))/(grid(2)-grid(1))

uBar

hatX = 0;

smoothIndShift = 3;

targetCell = 2;

linWgt = [3,1,1];

ru = MLWeno(basepolyncoeff,uBar,stencil,linWgt,targetCell,hatX,smoothIndShift,grid)

fun(0.3) - ru
