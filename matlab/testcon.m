% Test convergence of SWENO-AO in 1D
%
% Set up stencils (of cells) and linear weights, with the target cell 0.
% Set up the true solution with target point = 0.0.
% Input N so dx = 1/N is the element diameter.
% Reconstruct at target point 0.0, which is the left end of the target cell 0.
%
% Ex: Stencil [-2,1] is the mesh (-2,-1,0,1,2)*dx
% Specify jump value and place (give an integer, use times dx in the code)
function [] = testcon(dxInv,jumpValue,jumpPlace)
    if dxInv <= 0
        dxInv = 1;
    end
    dx = 1.0/dxInv;

    %% Set stencils
    %stencil = [ [-2,2]; [-2,0]; [-1,1]; [0,2] ];
    %linWgt  = [ 1/2 1/6 1/6 1/6];

    stencil = [ [-2,2]; [-2,0]; [-1,1];  [0,0] ];
    linWgt  = [ 5 2 3 1 ]; % normalize later

    %stencil = [ [-3,3]; [-2,0]; [-1,1];  [0,0] ];
    %linWgt  = [ 5 2 3 1 ]; % normalize later

    %% Set solution
    f = @(x) cos(x+0.2);
    
    % SETUP
    
    targetPt = 0.0;
    fJump = @(x) f(x) + jumpValue*(x<jumpPlace*dx);

    nStencil = size(stencil,1);
    
    for s = 1:nStencil
        if stencil(s,2) < stencil(s,1);
            st = stencil(s,1);
            stencil(s,1) = stencil(s,2);
            stencil(s,2) = st;
        end
        if stencil(s,1) > 0 || stencil(s,2) < 0
            warning('Each stencil must contain target cell 0.');
            return;
        end
    end
    
    nCells = max(stencil(:,2)) - min(stencil(:,1)) + 1;
    ic = 1 - min(stencil(:,1)); % index of target cell
    
    x = ( min(stencil(:,1)):1:max(stencil(:,2))+1 )*dx;

    % DEFINE UBAR AND INITIAL OUTPUT
    
    uBar = zeros(nCells,1);
    for i = 1:nCells
        uBar(i) = integral(f,x(i),x(i+1))/dx + jumpValue*(x(i)<jumpPlace*dx);
    end
    
    fprintf('UBar:');
    for i=1:length(uBar)
        fprintf(' %f',uBar(i));
    end
    fprintf('\n');
    
    fprintf('Stencils:              [');
    for s = 1:nStencil
        fprintf(' [%d %d]',stencil(s,1),stencil(s,2));
    end
    fprintf(' ]\n'); 

    % LINEAR WEIGHTS

    linWgt = linWgt / sum(linWgt);
    
    fprintf('Linear weights:       ');
    for i=1:length(linWgt)
        fprintf(' %f',linWgt(i));
    end
    fprintf('\n');    
    
    % SMOOTHNESS INDICATORS
    
    sigma = zeros(1,nStencil);
    for s = 1:nStencil
        if stencil(s,2) > stencil(s,1)
            for i = stencil(s,1):stencil(s,2)
                if i ~= 0
                    %sigma(s) = sigma(s) + ( ( uBar(i+ic) - uBar(ic) ) / i )^2;
                    sigma(s) = sigma(s) + ( uBar(i+ic) - uBar(ic) )^2;
                end
            end
            sigma(s) = sigma(s) / (stencil(s,2) - stencil(s,1));
            %midIndex = floor((stencil(s,1) + stencil(s,2))/2);
            %sigma(s) = sigma(s) + ( uBar(stencil(s,1)+ic) + uBar(stencil(s,2)+ic) ...
            %                        - 2*uBar(midIndex+ic) )^2;
        else
            sigma(s) = 0; %(dx/1)^2;
        end
    end

    fprintf('Smoothness indicators:');
    for s = 1:nStencil
        fprintf(' %f',sigma(s));
    end
    fprintf('\n');
    
    % NONLINEAR WEIGHTS

    hatWgt = zeros(1,nStencil);
    for s = 1:nStencil
        r = stencil(s,2) - stencil(s,1) + 1;
        eta = floor(r/2)+1;
        
        %hatWgt(s) = linWgt(s) / ( sigma(s)^eta + 1e-4*dx^r );
        %hatWgt(s) = linWgt(s) / ( sigma(s)^eta + (dx/5)^r );
        hatWgt(s) = linWgt(s) / ( ((r-1)*sigma(s)/4)^eta + (dx/5)^r );
    end
    nonlinWgt = hatWgt / sum(hatWgt);
    
    fprintf('Nonlinear weights:    ');
    for i=1:length(nonlinWgt)
        fprintf(' %f',nonlinWgt(i));
    end
    fprintf('\n');
    
    % THE RECONSTRUCTION

    targetVal = 0;
    for s = 1:nStencil
        basePolynCoeff = polyn(stencil(s,:));

        uBarStencil = uBar(stencil(s,1)+ic:stencil(s,2)+ic);
        targetVal = targetVal + nonlinWgt(s)*weno_reconst(targetPt/dx, basePolynCoeff, ...
                                                          uBarStencil, stencil(s,:)); 
    end
    
    % ERROR IN THE RECONSTRUCTION

    fprintf('Error = %5.3e\n',abs(targetVal - fJump(targetPt)));

end

% ==============================================================================

function [reconU] = weno_reconst(hatX, basePolynCoeff, uBarStencil, stencil)
    % hatX is relative coordinate x/dx

    r = stencil(2) - stencil(1) + 1;
   
    reconU = zeros(length(hatX));
    for k=1:r
        for p = 1:r
            reconU = reconU + uBarStencil(k) * basePolynCoeff(k,p)*hatX.^(p-1);
        end
    end
end

function [sol] = polyn(stencil)
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
