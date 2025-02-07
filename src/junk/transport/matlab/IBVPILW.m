function [] = IBVPILW(N,a,k)

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
%a = 1; k = 1;
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
    uBar = uBar + gaussWt(g)*init(a,k,gPt);
end

% Define reconstructions stencils
% A WENO (3,2) reconstruction
% Define reconstruction stencls
stencil32 = [[-1,1]; [-1,0]; [0,1]];
linWgt32  = [3,1,1];

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

% Define exact solution, initial and boundary conditions
% change it later for different conditions
%fexact = @(x,t) exp(-k*t)*sin(x-a*t);

%init = @(x) sin(x);

%boundaryL = @(t) exp(-k*t)*sin(-1-a*t);
%boundaryR = @(t) exp(-k*t)*sin( 1-a*t);

% Attach two small cells outside of the boundary
% In order to match with the physics boundary, following flow solver,
NTmax = 2000*a;
Tmax  = 0.5;
dt    = Tmax/NTmax;

CFL = a*dt/h;

NT    = NTmax;

uBarCurrent = uBar;
uBarNext = uBar;

diffRu = zeros(M,4);

alpha = 0.5; 
beta  = 1.5; 

% Inverse Lax-Wendroff scheme modefied for WENO finite 
% volume scheme, everything the same as NoGhostCell, expect
% replacing fixed boundary values with ILW approximation
figure, set(gcf)
%set(gca,'nextplot','replacechildren');
%filename = "ILW,a=" + num2str(a) +".avi";
%v = VideoWriter(filename);
%open(v);
for time = 1:NT
    currentT = time*dt;




    % exact solution
    %plot(linspace(0,1,100),fexact(linspace(0,1,100),currentT),'-')
    %plot(linspace(0,1,100),fexact(a,k,linspace(0,1,100),currentT),'-')
	 %hold on
    % Computational solution
    %plot(cell(1:end),uBarCurrent(1:end),'o');
    %[t,s] = title(['The Peclet number is ',num2str(Pe), ', CFL = ',num2str(CFL) ]);
	 %s.FontAngle = 'italic';
	 %legend({'Exact solution','Numerical solution'},'Location','northwest');
 
    %axis([0 1 -1 1])

    %frame = getframe(gcf);
	 %writeVideo(v,frame);

	 %pause(0.0001)

    %hold off

    %errorLnorm(2,uBarCurrent,a,k,x,currentT,1:N)


end

% ======================================================================
function [fu] = advectionFunc(a*u)

    % Linear advection case
    fu = a*u;

end

function [flux] = LaxFriedrich(a, uP, uM)

    flux = 0.5*(advectionFunc(uP) + advectionFunc(uM) - a*(uP-uM));

end

function [flux] = diffFlux(alpha,beta,ru)

    flux = ((ru(3)-ru(2))*beta^2/(2*alpha) - ...
            (ru(4)-ru(1))*alpha^2/(2*beta))/ ...
           (beta^2-alpha^2);

end

function [flux] = totalFlux(a, k, uP, uM, alpha, beta, ru, n)

    % Compute total flux consisting advection and diffusion flux
    % flux = au - kdu
    % n denotes the normal direction
    flux = (a*LaxFriedrich(a,uP,uM) - k*diffFlux(alpha,beta,ru))*n; 
end
