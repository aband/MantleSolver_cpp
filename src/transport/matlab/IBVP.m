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

% Weno (4,3) reconstruction for diffusive flux
stencil43 = [[-2,1];[-2,0];[-1,1]];
linWgt43  = [4,1,1]; 

stencil43L = [[0,3];[0,2];[0,1];[0,0]];
linWgt43L  = [4,1,1,1]; 

stencil43LL = [[-1,2];[-1,1];[0,2]];
linWgt43LL  = [4,1,1]; 

stencil43R = [[-3,0];[-2,0];[-1,0];[0,0]];
linWgt43R  = [4,1,1,1]; 

stencil43RR = [[-2,1];[-2,0];[-1,1]];
linWgt43RR  = [4,1,1]; 

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
    uBar = uBar + gaussWt(g)*sin(gPt);
end

% Biased WENO stencil method for Dirichlet boundary condition
eta_bias = [0,0,0,0];
%[uL, uR] = MultLWenoRecon(1.0, eta_bias, h,uBar,stencil32,linWgt32,2)

uBarCurrent = uBar;
uBarNext = uBar;

dt = Tmax/NT;

diffRu = zeros(M,4);

bVL = init(-1);
bVR = init(1);

alpha = 0.5; 
beta  = 1.5; 

% Time propogation and plotting
figure
for time = 0:Tmax:NT
    currentT = Tmax*time/NT;

    % Boundary treatment type 1
    UL = zeros(M,1);
    UM = zeros(M,1);

    [uL1, uR1] = MultLWenoRecon(1,eta_bias,h,uBarCurrent,stencil32L,linWgt32L,1);
    [uLN, uRN] = MultLWenoRecon(1,eta_bias,h,uBarCurrent,stencil32R,linWgt32R,N);

    UL(1) = init(-1);
    UR(1) = uL1; 
 
    UR(M) = init(1);
    UL(M) = uRN;

    UL(2) = uR1;
    UR(M-1) = uLN;

    for i=2:N-1
        [u1,u2] = MultLWenoRecon(1,eta_bias,h,uBarCurrent,stencil32,linWgt32,2);
        UR(i) = u1;
        UL(i+1) = u2;
    end

    % Calculate total flux and temperal update for uBar

    ru1 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43L,linWgt43L,alpha,beta,bVL,bVR,1);
    ru2 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43L,linWgt43L,alpha,beta,bVL,bVR,1);        
    uBarNext(1) = uBarCurrent(1) - dt/h*(totalFlux(a,UL(1),UR(1),alpha,beta,ru1,-1) +...
                                   totalFlux(a,UR(2),UL(2),alpha,beta,ru2, 1)) ;


    ru1 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43LL,linWgt43LL,alpha,beta,bVL,bVR,2);
    ru2 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43LL,linWgt43LL,alpha,beta,bVL,bVR,2);        
    uBarNext(2) = uBarCurrent(2) - dt/h*(totalFlux(a,UL(2),UR(2),alpha,beta,ru1,-1) +...
                                   totalFlux(a,UR(3),UL(3),alpha,beta,ru2, 1)) ;

    ru1 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43R,linWgt43R,alpha,beta,bVL,bVR,N);
    ru2 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43R,linWgt43R,alpha,beta,bVL,bVR,N);        
    uBarNext(N) = uBarCurrent(N) - dt/h*(totalFlux(a,UL(M-1),UR(M-1),alpha,beta,ru1,-1) +...
                                   totalFlux(a,UR(M),UL(M),alpha,beta,ru2, 1)) ;

    for s=3:N-1
        ru1 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43,linWgt43,alpha,beta,bVL,bVR,s);
        ru2 = MultLWenoReconVertx(1,eta_bias,h,uBarCurrent,stencil43,linWgt43,alpha,beta,bVL,bVR,s);        
        uBarNext(s) = uBarCurrent(s) - dt/h*(totalFlux(a,UL(s),UR(s),alpha,beta,ru1,-1) +...
                                       totalFlux(a,UR(s),UL(s),alpha,beta,ru2, 1)) ;
    end
   
    uBarCurrent = uBarNext;

    % exact solution
    plot(linspace(0,1,100),fexact(linspace(0,1,100),currentT),'-')
    plot(cell,uBarCurrent,'o');

    %axis([0 1 -1 1])
    pause(0.001)
end

[T,s] = title(['The Peclet number is ',num2str(Pe), ' ,h^{-1} is ',num2str(N)]);
s.FontAngle = 'italic';

% ================================================================================
function [fu] = advectionFunc(u)

    % Linear advection case
    fu = u;

end

function [flux] = LaxFriedrich(a, uP, uM)

    flux = 0.5*(advectionFunc(uP) + advectionFunc(uM) - a*(uP-uM));

end

function [flux] = diffFlux(alpha,beta,ru)

    flux = ((ru(3)-ru(2))*beta^2/alpha - ...
            (ru(4)-ru(1))*alpha^2/beta)/ ...
           (beta^2-alpha^2);

end

function [flux] = totalFlux(a, uP, uM, alpha, beta, ru, n)

    % Compute total flux consisting advection and diffusion flux
    % flux = au - kdu
    % n denotes the normal direction
    flux = (LaxFriedrich(a,uP,uM) - diffFlux(alpha,beta,ru))*n; 

end

% A new reconstruction function
function [uL,uR] = MultLWenoRecon(eps0, eta_bias, dx, uBar, stencil, linWgt, targetCell)

    nStencils = size(stencil,1);

    hatWgt = linWgt;

    maxR = max(stencil(:,2) - stencil(:,1) + 1);
    basePolynCoeff = zeros(nStencils,maxR,maxR);

    sigma = classicSmoothnessInd(uBar,stencil,targetCell);

    for s = 1:nStencils

        r = stencil(s,2) - stencil(s,1) + 1;

        left  = stencil(s,1) + targetCell;
        right = stencil(s,2) + targetCell;

        uBarStencil = uBar(left:right);

        eta = floor(r/2)+1;

        hatWgt(s) = linWgt(s) * ((sigma(s) + eps0*dx)/...
                                  (sigma(s) + (eps0*dx)^2))^r...
                               * (eps0*dx/(sigma(s)+eps0*dx))^eta_bias(s);

    end

    nonlinWgt = hatWgt / sum(hatWgt);

    uL = 0;
    uR = 0;

    for s=1:nStencils
        r = stencil(s,2) - stencil(s,1) + 1;

        basePolynCoeff(s,1:r,1:r) = polyn(stencil(s,:));
    end

    for s=1:nStencils

        r = stencil(s,2) - stencil(s,1) + 1;

        left  = stencil(s,1) + targetCell;
        right = stencil(s,2) + targetCell;

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

function [sol] = basePolynCoeff(stencil)
    r = stencil(2) - stencil(1) + 1;

    M = zeros(r,r);
    for j = 1:r
        xLeft  = stencil(1) + j - 1;
        xRight = stencil(1) + j;
    
        for p=1:r
            M(j,p) = xRight^p/p - xLeft^p/p;
        end
    end

    sol = zeros(r,r);
    
    for k = 1:r
        B = zeros(r,1);
        B(k) = 1;
        sol(k,:) = M\B;
    end
end

function [val] = polynEval(hatX, basePolynCoeff, uBarStencil, stencil)
    % hatX is relative coordinate x/dx

    r = stencil(2) - stencil(1) + 1;
   
    val = zeros(length(hatX));
    for k=1:r
        for p = 1:r
            val = val + uBarStencil(k) * basePolynCoeff(k,p)*hatX.^(p-1);
        end
    end
end

function [val] = polynEvalDer(ell, hatX, basePolynCoeff, uBarStencil, stencil)
    % hatX is relative coordinate x/dx

    r = stencil(2) - stencil(1) + 1;
   
    val = zeros(1,length(hatX));
    for k=1:r
        for p = 1+ell:r
            val = val + uBarStencil(k) ...
                  * basePolynCoeff(k,p)*hatX.^(p-1-ell)*factorial(p-1)/factorial(p-1-ell);
        end
    end
end

function [sigma] = classicSmoothnessInd(uBar, stencil, ic)
    gaussPt = [ -sqrt(5 + 2*sqrt(10/7))/3 -sqrt(5 - 2*sqrt(10/7))/3 ...
                0 sqrt(5 - 2*sqrt(10/7))/3 sqrt(5 + 2*sqrt(10/7))/3 ];
    gaussWt = [ (322 - 13*sqrt(70))/900 (322 + 13*sqrt(70))/900 ...
                128/225 (322 + 13*sqrt(70))/900 (322 - 13*sqrt(70))/900 ];
    
    gaussPt = (gaussPt + 1)/2;
    gaussWt = gaussWt/2;

    nStencil = size(stencil,1);
    sigma = zeros(nStencil,1);
    
    for s=1:nStencil
        r = stencil(s,2) - stencil(s,1) + 1;
        uBarStencil = uBar(stencil(s,1)+ic:stencil(s,2)+ic);
        coeff = basePolynCoeff(stencil(s,:));

        for ell=1:r-1 %min(2,r-1)
            val = polynEvalDer(ell, gaussPt, coeff, uBarStencil, stencil(s,:));
            sigma(s) = sigma(s) + (val.*val)*gaussWt';
        end    
    end
end

% =============================================================================================

function [ru] = MultLWenoReconVertx(eps0, eta_bias, dx, uBar, stencil, linWgt, alpha, beta, bVL, bVR, targetCell)

    N = size(uBar,1);

    nStencils = size(stencil,1);

    hatWgt = linWgt;

    maxR = max(stencil(:,2) - stencil(:,1) + 1);
    basePolynCoeff = zeros(nStencils,maxR,maxR);

    sigma = classicSmoothnessInd(uBar,stencil,targetCell);

    for s = 1:nStencils

        r = stencil(s,2) - stencil(s,1) + 1;

        left  = stencil(s,1) + targetCell;
        right = stencil(s,2) + targetCell;

        uBarStencil = uBar(left:right);

        eta = floor(r/2)+1;

        hatWgt(s) = linWgt(s) * ((sigma(s) + eps0*dx)/...
                                 (sigma(s) + (eps0*dx)^2))^r...
                               * (eps0*dx/(sigma(s)+eps0*dx))^eta_bias(s);

    end

    nonlinWgt = hatWgt / sum(hatWgt);

    ru = zeros(4,1);

    for s=1:nStencils
        r = stencil(s,2) - stencil(s,1) + 1;

        basePolynCoeff(s,1:r,1:r) = polyn(stencil(s,:));
    end

    for s=1:nStencils

        r = stencil(s,2) - stencil(s,1) + 1;

        left  = stencil(s,1) + targetCell;
        right = stencil(s,2) + targetCell;

        uBarStencil = uBar(left:right);

        % Create coefficients of polynomial
        r = right - left + 1;

		  if targetCell == 1

				ru(1) = bVL;
				ru(2) = bVL;
				ru(3) = ru(3) + nonlinWgt(s)*weno_reconst(alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(4) = ru(4) + nonlinWgt(s)*weno_reconst(beta-0.5,basePolynCoeff,s,uBarStencil,r);

		  elseif targetCell == N

				ru(1) = ru(1) + nonlinWgt(s)*weno_reconst(-beta-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(2) = ru(2) + nonlinWgt(s)*weno_reconst(-alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(3) = bVR;
				ru(4) = bVR;

		  elseif targetCell == 2

				ru(1) = bVL;
				ru(2) = ru(1) + nonlinWgt(s)*weno_reconst(-alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(3) = ru(2) + nonlinWgt(s)*weno_reconst(alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(4) = ru(3) + nonlinWgt(s)*weno_reconst(beta-0.5,basePolynCoeff,s,uBarStencil,r);

        elseif targetCell == N-1

				ru(1) = ru(1) + nonlinWgt(s)*weno_reconst(-beta-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(2) = ru(2) + nonlinWgt(s)*weno_reconst(-alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(3) = ru(3) + nonlinWgt(s)*weno_reconst(alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(4) = bVR;

        else

				ru(1) = ru(1) + nonlinWgt(s)*weno_reconst(-beta-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(2) = ru(2) + nonlinWgt(s)*weno_reconst(-alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(3) = ru(3) + nonlinWgt(s)*weno_reconst(alpha-0.5,basePolynCoeff,s,uBarStencil,r);
				ru(4) = ru(4) + nonlinWgt(s)*weno_reconst(beta-0.5,basePolynCoeff,s,uBarStencil,r);

        end

    end

end
