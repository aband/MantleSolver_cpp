% 1D code utilizing WENO reconstruction for FV scheme
% A Initial-Boundary value problem

clc; clear

% Define cell number
N = 64;

% Define vertex number
M = N + 1;

% Define right hand side
x = linspace(-1,1,M);
h = 1/N;
cell = linspace(-1,1,N);
sol = zeros(size(cell));

%f = 4*(x(2:end).^3 - x(1:end-1).^3);

% Define convection and diffusion coefficients
% A convection dominated convection-diffusion problem
a = 200; k = 1;
Pe = a/k;

% A WENO (3,2) reconstruction
% Define reconstruction stencls
stencil32 = [[-1,1]; [-1,0]; [0,1]];
linWgt32  = [3,1,1];

% Degenerate multilevel WENO (3,2) - (2,1) reconstruction at boundary
stencil21L = [[0,1];[0,0]];
linWgt21L  = [2,1];

stencil21R = [[-1,0];[0,0]];
linWgt21R  = [2,1];

% Non degenerate biased multilevel WENO (3,2) reconstruction at boundary
% Mimicing biased finite difference approximation
% Multilevel WENO can be used to deal with discontinuity appearing near the boundary 
% by going down to order zero!
stencil32L = [[0,2];[0,1];[0,0]];
linWgt32L  = [3,2,1];

stencil32R = [[-2,0];[-1,0];[0,0]];
linWgt32R  = [3,2,1];

% Define exact solution, initial and boundary conditions
% change it later for different conditions
fexact = @(x,t) exp(-k*t)*sin(x-a*t);

init = @(x) sin(x);

boundaryL = @(t) exp(-k*t)*sin(-1-a*t);
boundaryR = @(t) exp(-k*t)*sin( 1-a*t);

% Attach two small cells outside of the boundary
% In order to match with the physics boundary, following flow solver,
NT   = 1000;
Tmax = 0.5;
dt   = Tmax/NT;

% Lax-Friedrich type numerical flux
LF = a;
hatF = @(a,b) 0.5*( f(a) + f(b) - LF*(b-a) );

% Gauss quadrature points and weights
gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
gaussWt = [5/18 8/18 5/18]; %half value!

% Create uBar
uBar = 0;

vertxL = x(1:end-1); 
vertxR = x(2:end);
for g = 1:3
    gPt = vertxL + (vertxR-vertxL)*(gaussPt(g)+1)/2;
    uBar = uBar + gPt;
end

% Biased WENO stencil method for Dirichlet boundary condition
[uL, uR] = MultLWenoRecon(1.0,h,uBar,stencil32,linWgt32,2);
init(cell(2))


% Time propogation and plotting
figure
for time = 0:Tmax:NT
    currentT = Tmax*time/NT;

    
    % exact solution
    plot(linspace(0,1,100),fexact(linspace(0,1,100),currentT),'-')
    axis([0 1 -1 1])
    pause(0.01)
end
[T,s] = title(['The Peclet number is ',num2str(Pe), ' ,h^{-1} is ',num2str(N)]);
s.FontAngle = 'italic';

% ================================================================================
% A new reconstruction function
function [uL,uR] = MultLWenoRecon(eps0, dx, uBar, stencil, linWgt, targetCell)

    nStencils = size(stencil,1);

    hatWgt = linWgt;

    maxR = max(stencil(:,2) - stencil(:,1) + 1);
    basePolynCoeff = zeros(nStencils,maxR,maxR);

    for s = 1:nStencils

        r = stencil(s,2) - stencil(s,1) + 1;

        left  = stencil(s,1) + targetCell;
        right = stencil(s,2) + targetCell;

        uBarStencil = uBar(left:right);

        sigma = computeSigma(uBarStencil,stencil(s,:));

        eta = floor(r/2)+1;

        hatWgts(s) = linWgt(s) / (sigma^eta + eps0*(dx/100)^r);

    end

    nonlinWgt = hatWgt / sum(hatWgt);

    uL = 0;
    uR = 0;

    for s=1:nStencils

        r = stencil(s,2) - stencil(s,1) + 1;

        basePolynCoeff(s,1:r,1:r) = polyn(stencil(s,:));

        left  = stencil(s,1) + targetCell;
        right = stencil(s,1) + targetCell;

        uBarStencil = uBar(left:right);

        % Create coefficients of polynomial
        r = right - left + 1;

        uL = uL + nonlinWgt(s)*weno_reconst(-0.5,basePolynCoeff,s,uBarStencil,r); 
        uR = uR + nonlinWgt(s)*weno_reconst( 0.5,basePolynCoeff,s,uBarStencil,r); 
    end

end

function [ru] = weno_reconst(p,coeff,s,uBarStencil,r)

   % p is relative coordinate (x-x0)/dx

   ru = 0.0;

   for k=1:r
       for i = 1:r
           ru = ru + uBarStencil(k) * coeff(s,k,i)*p.^(i-1);
       end
   end

end

function [sol] = polyn(stencil)

r = stencil(2) - stencil(1) + 1;

M = zeros(r,r);

sol = zeros(r,r);

    for k = 1:r
       for j = 1:r
           xleft  = stencil(1)+(j-1)-0.5;
           xright = xleft+1;

           for i=1:r
               M(j,i) = integral(@(x) x.^(i-1), xleft, xright);
           end
       end

       B = zeros(r,1);
       B(k) = 1;
       sol(k,:) = M\B;
    end

end

function [Sigma] = computeSigma(uBar, stencil)

    Sigma = 0;

    target = 1-stencil(1);

    if stencil(2)> stencil(1)
        for i = stencil(1):stencil(2)
            if i~= 0
                I = i-stencil(1) + 1;
                Sigma = Sigma + (uBar(I) - uBar(target))^2;
                %Sigma = Sigma + ( ( uBar(I) - uBar(target) ) / i )^2;
            end
        end
        %Sigma = Sigma /(stencil(2) - stencil(1));
        %midIndex = floor((stencil(1) + stencil(2))/2) - stencil(1) + 1;
        %Sigma = Sigma + ( uBar(1) + uBar(end) - 2*uBar(midIndex) )^2;
 
    else
        Sigma = 0;
    end

end

function [Sigma] = computeSigma2(uBar, stencil)

    Sigma = 0;

    target = 1-stencil(1);

    if stencil(2)> stencil(1)
        for i = stencil(1):stencil(2)
            %if i~= 0
						  for j = stencil(1):stencil(2)
                I = i-stencil(1) + 1;
					 J = j-stencil(1) + 1;
                %Sigma = Sigma + (uBar(I) - uBar(target)^2);
                %Sigma = Sigma + ( ( uBar(I) - uBar(target) ) / i )^2;
					 if J ~= I
				    Sigma = Sigma + ( (uBar(I) - uBar(J) ) / (J-I) )^2;
					 end
						  end
            %end
        end
		  Sigma = Sigma / ((stencil(2)-stencil(1))*(stencil(2)-stencil(1)+1)/2);
        %Sigma = Sigma /(stencil(2) - stencil(1));
        %midIndex = floor((stencil(1) + stencil(2))/2) - stencil(1) + 1;
        %Sigma = Sigma + ( uBar(1) + uBar(end) - 2*uBar(midIndex) )^2;
 
    else
        Sigma = 0;
    end

end



