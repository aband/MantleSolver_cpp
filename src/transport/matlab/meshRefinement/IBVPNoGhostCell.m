function [] = IBVPNoGhostCell(N,a,k)

% Define vertex number
M = N + 1;

% Define right hand side
startP = -1;
endP = 1;


refinement = 20;
href = 0.2/refinement;
x = [linspace(-1,-0.8,refinement),linspace(-0.8+href,0.8-href,M-2*refinement),linspace(0.8,1,refinement)];

%x = linspace(startP,endP,M)
%size(x)

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
NTmax = 100*a;
Tmax  = 0.5;
dt    = Tmax/NTmax;

CFL = a*dt/h;

NT    = NTmax;

uBarCurrent = uBar;
uBarNext = uBar;

diffRu = zeros(M,4);

alpha = 0.5; 
beta  = 1.5; 

% Calculate base polynomial coefficients in advance
nStencils32 = size(stencil32,1);
nStencils43 = size(stencil43,1); 

maxR32 = max(stencil32(:,2)-stencil32(:,1)+1);
maxR43 = max(stencil43(:,2)-stencil43(:,1)+1);

basepolyncoeff32All = zeros(N-2,nStencils32,maxR32,maxR32);
basepolyncoeff43All = zeros(N-3,nStencils43,maxR43,maxR43);

basepolyncoeff32All(1,:,:,:) = basePolynRefine(2,stencil32,-0.5,x);

for s = 3:N-1
    basepolyncoeff32All(s-1,:,:,:) = basePolynRefine(s,stencil32,-0.5,x);
    basepolyncoeff43All(s-2,:,:,:) = basePolynRefine(s,stencil43, 0.0,x);
end

basepolyncoeff32L = basePolynRefine(1,stencil32L,-0.5,x);
basepolyncoeff43L = basePolynRefine(1,stencil43L, 0.0,x);
basepolyncoeff43LL = basePolynRefine(2,stencil43LL, 0.0,x);

basepolyncoeff32R = basePolynRefine(N,stencil32R,-0.5,x);
basepolyncoeff43R = basePolynRefine(N,stencil43R, 0.0,x);
basepolyncoeff43RR = basePolynRefine(N,stencil43R, -1.0,x);

% Time propogation and plotting
% Forward Eurlar time propogation

clf;
drawnow;
for time = 1:NT
    currentT = time*dt;

    %bVL = boundaryL(currentT);
    %bVR = boundaryR(currentT);

    [bVL,bVR] = boundary(a,k,currentT);

    % Update interior cells first 
    for s = 3:N-2
        uLp = multiLWENORefine(reshape(basepolyncoeff32All(s-2,:,:,:),nStencils32,maxR32,maxR32),h,uBarCurrent,stencil32,linWgt32,s-1, 0.5,3,x); 
        uLm = multiLWENORefine(reshape(basepolyncoeff32All(s-1,:,:,:),nStencils32,maxR32,maxR32),h,uBarCurrent,stencil32,linWgt32,s  ,-0.5,3,x);
        uRm = multiLWENORefine(reshape(basepolyncoeff32All(s-1,:,:,:),nStencils32,maxR32,maxR32),h,uBarCurrent,stencil32,linWgt32,s  , 0.5,3,x);
        uRp = multiLWENORefine(reshape(basepolyncoeff32All(s,:,:,:),nStencils32,maxR32,maxR32),h,uBarCurrent,stencil32,linWgt32,s+1,-0.5,3,x); 

        hatX = [-1*beta,-1*alpha,alpha,beta];

        ruL = multiLWENORefine(reshape(basepolyncoeff43All(s-2,:,:,:),nStencils43,maxR43,maxR43),h,uBarCurrent,stencil43,linWgt43,s  ,hatX,2,x);
        ruR = multiLWENORefine(reshape(basepolyncoeff43All(s-1,:,:,:),nStencils43,maxR43,maxR43),h,uBarCurrent,stencil43,linWgt43,s+1,hatX,2,x);

        uBarNext(s) = uBarCurrent(s) - dt/h * (totalFlux(a,k,uLm,uLp,alpha*h,beta*h,ruL,-1) +...
                                               totalFlux(a,k,uRp,uRm,alpha*h,beta*h,ruR, 1));

    end

    % Treat boundary without assigning ghost cells
    % Left boundary cell
    s = 1;     
    uLp = bVL;
    uLm = multiLWENORefine(basepolyncoeff32L,h,uBarCurrent,stencil32L,linWgt32L,1,-0.5,3,x);
    uRm = multiLWENORefine(basepolyncoeff32L,h,uBarCurrent,stencil32L,linWgt32L,1, 0.5,3,x);
    uRp = multiLWENORefine(reshape(basepolyncoeff32All(1,:,:,:),nStencils32,maxR32,maxR32), h,uBarCurrent,stencil32 ,linWgt32 ,2,-0.5,1,x);

    hatX = [alpha,beta];
    ruL = multiLWENORefine(basepolyncoeff43L,h,uBarCurrent,stencil43L ,linWgt43L ,s,hatX,1,x);
    ruL = [bVL,bVL,ruL];

    %hatX = [-1*alpha,alpha,beta];
    %ruR = multiLWENORefine(reshape(basepolyncoeff43LL,h,uBarCurrent,stencil43LL,linWgt43LL,s+1,hatX,2);
    %ruR = [bVL,ruR];

    hatX = [-1*beta,-1*alpha,alpha,beta]/beta;
    ruR = multiLWENORefine(basepolyncoeff43LL,h,uBarCurrent,stencil43LL,linWgt43LL,s+1,hatX,2,x);

    %leftflux = totalFlux(a,k,uLp,uLm,alpha*h,beta*h,ruL,-1);
    leftflux = -80.0; 

	 rightflux = totalFlux(a,k,uRp,uRm,alpha*h,beta*h,ruR, 1);

    uBarNext(s) = uBarCurrent(s) - dt/h * (leftflux +...
                                           rightflux);

    s = 2; % The second cell
    uLp = multiLWENORefine(basepolyncoeff32L,h,uBarCurrent,stencil32L,linWgt32L,s-1, 0.5,3,x); 
    uLm = multiLWENORefine(reshape(basepolyncoeff32All(1,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s  ,-0.5,1,x);
    uRm = multiLWENORefine(reshape(basepolyncoeff32All(1,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s  , 0.5,1,x);
    uRp = multiLWENORefine(reshape(basepolyncoeff32All(2,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s+1,-0.5,1,x);

    %hatX = [-1*alpha,alpha,beta];
    %ruL = multiLWENORefine(reshape(basepolyncoeff43LL,h,uBarCurrent,stencil43LL,linWgt43LL,s,hatX,2);
    %ruL = [bVL,ruL];

    hatX = [-1*beta,-1*alpha,alpha,beta]/beta;
    ruL = multiLWENORefine(basepolyncoeff43LL,h,uBarCurrent,stencil43LL,linWgt43LL,s,hatX,2,x);

    hatX = [-1*beta,-1*alpha,alpha,beta];
    ruR = multiLWENORefine(reshape(basepolyncoeff43All(s-1,:,:,:),nStencils43,maxR43,maxR43)  ,h,uBarCurrent,stencil43  ,linWgt43  ,s+1,hatX,2,x);

    uBarNext(s) = uBarCurrent(s) - dt/h * (totalFlux(a,k,uLm,uLp,alpha*h,beta*h,ruL,-1) +...
                                           totalFlux(a,k,uRp,uRm,alpha*h,beta*h,ruR, 1));

    % Right boundary cell  
    s = N;
    uLp = multiLWENORefine(reshape(basepolyncoeff32All(s-2,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s-1, 0.5,3,x); 
    uLm = multiLWENORefine(basepolyncoeff32R,h,uBarCurrent,stencil32R,linWgt32R,s  ,-0.5,3,x);
    uRm = multiLWENORefine(basepolyncoeff32R,h,uBarCurrent,stencil32R,linWgt32R,s  , 0.5,3,x);
    uRp = bVR; 

    %hatX = [-1*beta,-1*alpha,alpha];
    %ruL = multiLWENORefine(reshape(basepolyncoeff43R,h,uBarCurrent,stencil43R,linWgt43R,s,hatX,2);
    %ruL = [ruL,bVR];

    hatX = [-1*beta,-1*alpha,alpha,beta]/beta;
    ruL = multiLWENORefine(basepolyncoeff43R,h,uBarCurrent,stencil43R,linWgt43R,s,hatX,2,x);

    hatX = [-1*beta,-1*alpha];
    ruR = multiLWENORefine(basepolyncoeff43RR,h,uBarCurrent,stencil43R,linWgt43R,s,hatX,1,x);
    ruR = [ruR,bVR,bVR];

    leftfluxN = totalFlux(a,k,uLm,uLp,alpha*h,beta*h,ruL,-1);
    
	 rightfluxN = totalFlux(a,k,uRp,uRm,alpha*h,beta*h,ruR, 1);
    %rightfluxN = 15;

    uBarNext(s) = uBarCurrent(s) - dt/h * (leftfluxN +...
                                           rightfluxN);

    % The second last cell
    s = N-1;
    uLp = multiLWENORefine(reshape(basepolyncoeff32All(s-2,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s-1, 0.5,3,x); 
    uLm = multiLWENORefine(reshape(basepolyncoeff32All(s-1,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s  ,-0.5,3,x);
    uRm = multiLWENORefine(reshape(basepolyncoeff32All(s-1,:,:,:),nStencils32,maxR32,maxR32) ,h,uBarCurrent,stencil32 ,linWgt32 ,s  , 0.5,3,x);
    uRp = multiLWENORefine(basepolyncoeff32R,h,uBarCurrent,stencil32R,linWgt32R,s+1,-0.5,3,x); 

    hatX = [-1*beta,-1*alpha,alpha,beta];
    ruL = multiLWENORefine(reshape(basepolyncoeff43All(s-2,:,:,:),nStencils43,maxR43,maxR43),h,uBarCurrent,stencil43,linWgt43,s,hatX,2,x);

    %hatX = [-1*beta,-1*alpha,alpha];
    %ruR = multiLWENORefine(reshape(basepolyncoeff43R,h,uBarCurrent,stencil43R,linWgt43R,s+1,hatX,2);
    %ruR = [ruR,bVR];

    hatX = [-1*beta,-1*alpha,alpha,beta]/beta;
    ruR = multiLWENORefine(basepolyncoeff43R,h,uBarCurrent,stencil43R,linWgt43R,s+1,hatX,2,x);

    uBarNext(s) = uBarCurrent(s) - dt/h * (totalFlux(a,k,uLm,uLp,alpha*h,beta*h,ruL,-1) +...
                                           totalFlux(a,k,uRp,uRm,alpha*h,beta*h,ruR, 1));

    % Update current uBar with next uBar
    uBarCurrent = uBarNext;

    % exact solution
	 %if mod(time,100) == 0
        clf;
        %plot(linspace(-1,1,100),fexact(a,k,linspace(-1,1,100),currentT),'-')
	     %hold on
        plot(cell(1:end),uBarCurrent(1:end),'o');
        [t,s] = title(['Pe =  ',num2str(Pe), ', CFL = ',num2str(CFL) ,', Time = ', num2str(currentT)...
                       ', N = ',num2str(N)]);
	     s.FontAngle = 'italic';
	     %legend({'Exact solution','Numerical solution'},'Location','northwest');
 
        axis([-1 1 -0.1 1.1])

        %frame = getframe(gcf);
	     %writeVideo(v,frame);

        %errorLnorm(2,uBarCurrent,a,k,x,currentT,1:N)

        %pause(0.0001)
	     drawnow;
		  %pause;
        %hold off
    %end

end
%close(v);


end

% ======================================================================
function [fu] = advectionFunc(a,u)

    % Linear advection case
    fu = a*u;

end

function [flux] = LaxFriedrich(a, uP, uM)

    flux = 0.5*(advectionFunc(a,uP) + advectionFunc(a,uM) - abs(a)*(uP-uM));

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
    flux = (1*LaxFriedrich(a,uP,uM) - k*diffFlux(alpha,beta,ru))*n; 
end
