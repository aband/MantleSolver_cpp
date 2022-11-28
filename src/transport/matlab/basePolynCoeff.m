function [basepolyncoeff] = basePolynCoeff(stencil,centerShift)

nStencils = size(stencil,1);

maxR = max(stencil(:,2)-stencil(:,1)+1);
basepolyncoeff = zeros(nStencils,maxR,maxR);

r = zeros(nStencils,1);

for s = 1:nStencils
   r(s) = stencil(s,2) - stencil(s,1) + 1;
   basepolyncoeff(s,1:r(s),1:r(s)) = polyn(stencil(s,:),centerShift);
end

end

function [sol] = polyn(stencil,center)
    r = stencil(2) - stencil(1) + 1;

    M = zeros(r,r);
    for j = 1:r
        xLeft  = stencil(1) + j - 1 + center;
        xRight = stencil(1) + j + center;
    
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


