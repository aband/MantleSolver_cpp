function [basepolyncoeff] = basePolynRefine(targetCell,stencil,centerShift,x)

nStencils = size(stencil,1);

maxR = max(stencil(:,2)-stencil(:,1)+1);
basepolyncoeff = zeros(nStencils,maxR,maxR);

r = zeros(nStencils,1);

for s = 1:nStencils
   r(s) = stencil(s,2) - stencil(s,1) + 1;
   basepolyncoeff(s,1:r(s),1:r(s)) = polyn(targetCell,stencil(s,:),centerShift,x);
end

end

function [sol] = polyn(targetCell,stencil,center,x)
    r = stencil(2) - stencil(1) + 1;

    M = zeros(r,r);

    % Extract target cell length
    cLeftInd = targetCell;
    cRightInd = targetCell+1;

    hCenter = x(cRightInd)-x(cLeftInd);

    cCenter = x(cLeftInd) - hCenter*center;

    for j = 1:r
        % Get left and right vertx index
        vLeftInd = targetCell + stencil(1) + j-1;
        vRightInd = targetCell + stencil(1) + j;

        xLeft = x(vLeftInd);
        xRight = x(vRightInd);

        %xLeftold  = stencil(1) + j - 1 + center
        %xRightold = stencil(1) + j + center
   
        xLeftShift = (xLeft-cCenter)/hCenter;
        xRightShift = (xRight-cCenter)/hCenter;

        for p=1:r
            M(j,p) = xRightShift^p/p - xLeftShift^p/p;
        end
    end

    sol = zeros(r,r);

    hstencil = targetCell+(stencil(1):stencil(2));

    hHat = (x(hstencil+1)-x(hstencil))/hCenter;

    for k = 1:r
        B = zeros(r,1);
        B(k) = hHat(k);
        sol(k,:) = M\B;
    end
end
