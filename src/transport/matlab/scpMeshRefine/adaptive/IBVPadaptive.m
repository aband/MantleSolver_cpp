function [] = IBVPadaptive(N,a,k,CFL,threshold)

% Define right hand side
startP = -1;
endP = 1;

% Define base grid x
%baseGrid = linspace(startP,endP,M)

baseGrid = generateMesh(threshold,N);

M = length(baseGrid);
N = M-1;

cell = (baseGrid(1:end-1)+baseGrid(2:end))/2;
sol = zeros(size(cell));

baseh = (endP-startP)/N;

% Define convection and diffusion coefficients
% A convection dominated convection-diffusion problem
%a = 1; k = 1;
Pe = a/k;

% Gauss quadrature points and weights
gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
gaussWt = [5/18 8/18 5/18]; %half value!

% Create uBar
uBar = 0;

vertxL = baseGrid(1:end-1); 
vertxR = baseGrid(2:end);
for g = 1:3
    gPt = vertxL + (vertxR-vertxL)*(gaussPt(g)+1)/2;
    uBar = uBar + gaussWt(g)*init(a,k,gPt);
end

% Define time integral parameters
Tmax  = 0.5;
dt    = CFL*baseh/(2*a);
NTmax = ceil(Tmax/dt);
dt    = Tmax/NTmax;

NT    = NTmax;

uBarCurrent = uBar;
uBarNext = uBar;

diffRu = zeros(M,4);

alpha = 0.5; 
beta  = 1.5; 

refineLevel = 2;

% Define stencils
global stencil32
global linWgt32  
global stencil32L
global linWgt32L 
global stencil32R 
global linWgt32R  
global stencil43
global linWgt43  
global stencil43L 
global linWgt43L  
global stencil43LL
global linWgt43LL  
global stencil43R
global linWgt43R   
global nStencils32
global nStencils43 
global maxR32 
global maxR43 

stencil32 = [[-1,1]; [-1,0]; [0,1]];
linWgt32  = [2,1,1];

stencil32L = [[0,2];[0,1];[0,0]];
linWgt32L  = [2,1,1];

stencil32R = [[-2,0];[-1,0];[0,0]];
linWgt32R  = [2,1,1];

% Weno (4,3) reconstruction for diffusive flux
stencil43 = [[-2,1];[-2,0];[-1,1];[0,0]];
linWgt43  = [4,1,1,1]; 

stencil43L = [[0,3];[0,2];[0,1];[0,0]];
linWgt43L  = [4,1,1,1]; 

stencil43LL = [[-1,2];[-1,1];[0,2];[0,0]];
linWgt43LL  = [4,1,1,1]; 

stencil43R = [[-3,0];[-2,0];[-1,0];[0,0]];
linWgt43R  = [4,1,1,1]; 


% Create stencil polynomials with given grid information

nStencils32 = size(stencil32,1);
nStencils43 = size(stencil43,1); 

maxR32 = max(stencil32(:,2)-stencil32(:,1)+1);
maxR43 = max(stencil43(:,2)-stencil43(:,1)+1);

uBarCurrent = uBar;
uBarNext = uBar;

[stencilPolyn32,stencilPolyn43] = createStencilPolyn(baseGrid);

% Time propogation starting here
% Forward Eurlar time propogation

refineLevel = 2;

clf;
drawnow;
for time = 1:NT
    currentT = time*dt;

    [bVL,bVR] = boundary(a,k,currentT);

    % Create refined grid, refined uBar and corresponding stencil polynomials
	 % For each time step
    [childGrid, childUBar, loc] = refineMesh(baseGrid, uBarCurrent, threshold, refineLevel, stencilPolyn32, stencilPolyn43);

    gridNow = childGrid;
    uBarCurrent = childUBar;
    [stencilPolyn32,stencilPolyn43] = createStencilPolyn(gridNow);

    uBarNext = uBarCurrent;

    Nnow = length(uBarCurrent)

    % Create stencil polynomials
    uR = [bVL];
    uL = [];
    for s = 1:Nnow

        if s == 1
           stencil = stencil32L;
           linWgt  = linWgt32L;
        elseif s == Nnow
           stencil = stencil32R;
           linWgt  = linWgt32R;
        else
           stencil = stencil32;
           linWgt  = linWgt32;
        end 

        ur = multiLWENORefine(reshape(stencilPolyn32(s,:,:,:),nStencils32,maxR32,maxR32),uBarCurrent,stencil,linWgt,s, 0.5,3,gridNow); 
        ul = multiLWENORefine(reshape(stencilPolyn32(s,:,:,:),nStencils32,maxR32,maxR32),uBarCurrent,stencil,linWgt,s  ,-0.5,3,gridNow);

        uR = [uR,ur];
        uL = [uL,ul];

	 end

    uL = [uL,bVR]; 

    for s = 1:Nnow
 
        ruL = getru(bVL,bVR,alpha,beta,s  ,reshape(stencilPolyn43(s,:,:,:)  ,nStencils43,maxR43,maxR43),gridNow,uBarCurrent);
        ruR = getru(bVL,bVR,alpha,beta,s+1,reshape(stencilPolyn43(s+1,:,:,:),nStencils43,maxR43,maxR43),gridNow,uBarCurrent);

        uLp = uR(s);
        uLm = uL(s);

        uRm = uR(s+1);
        uRp = uL(s+1);

        h = gridNow(s+1) - gridNow(s);

        uBarNext(s) = uBarCurrent(s) - dt/h * (totalFlux(a,k,uLm,uLp,alpha*h,beta*h,ruL,-1) +...
                                               totalFlux(a,k,uRp,uRm,alpha*h,beta*h,ruR, 1));

    end


    % Update current uBar with next uBar
    length(uBarNext)
	 length(uBarCurrent)
    length(gridNow)

    uBarCurrent = uBarNext;
  


    % exact solution
	 %if mod(time,100) == 0
        clf;
        %plot(linspace(-1,1,100),fexact(a,k,linspace(-1,1,100),currentT),'-')
	     %hold on
		  cellNow = (gridNow(1:end-1)+gridNow(2:end))/2;
        plot(cellNow,uBarCurrent,'o');
		  %hold on
		  %refinedCell = (childGrid(1:end-1)+childGrid(2:end))/2;
		  %plot(refinedCell,childUBar,'*');
        %plot(cell,coarseUBar,'-+');
		  %hold off
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
		  pause;
        %hold off
    %end

    childUBar = uBarCurrent;

    [coarseUBar] = coarseMesh(childGrid, childUBar, loc, refineLevel);

    uBarCurrent = coarseUBar;

end
%close(v);


end

% ======================================================================
function [ru] = getru(bVL,bVR,alpha,beta,targetCell, basepolyncoeff,baseGrid,uBarCurrent)

N = length(uBarCurrent);

global stencil32
global linWgt32  
global stencil32L
global linWgt32L 
global stencil32R 
global linWgt32R  
global stencil43
global linWgt43  
global stencil43L 
global linWgt43L  
global stencil43LL
global linWgt43LL  
global stencil43R
global linWgt43R   
global nStencils32
global nStencils43 
global maxR32 
global maxR43 

if targetCell == 1
    hatX = [alpha,beta];
    ru   = multiLWENORefine(basepolyncoeff,uBarCurrent,stencil43L ,linWgt43L ,targetCell,hatX,1,baseGrid);
    ru   = [bVL,bVL,ru];
elseif targetCell == 2
    hatX = [-1*beta,-1*alpha,alpha,beta]/beta;
    ru   = multiLWENORefine(basepolyncoeff,uBarCurrent,stencil43LL,linWgt43LL,targetCell,hatX,2,baseGrid);
elseif targetCell == N
    hatX = [-1*beta,-1*alpha,alpha,beta]/beta;
    ru   = multiLWENORefine(basepolyncoeff,uBarCurrent,stencil43R,linWgt43R,targetCell,hatX,2,baseGrid);
elseif targetCell == N+1
    hatX = [-1*beta,-1*alpha];
    ru   = multiLWENORefine(basepolyncoeff,uBarCurrent,stencil43R,linWgt43R,targetCell-1,hatX,1,baseGrid);
    ru   = [ru,bVR,bVR];
else
    hatX = [-1*beta,-1*alpha,alpha,beta];
    ru = multiLWENORefine(basepolyncoeff,uBarCurrent,stencil43,linWgt43,targetCell,hatX,2,baseGrid);

end

end


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

