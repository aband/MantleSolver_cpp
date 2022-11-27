% 1D Point-wise multi level weno reconstruction

function [ru] = multiLWENO1D(x,dx,uBar,stencil,linWgt,targetCell,hatX,centerShift,smoothIndShift)

nStencils = size(stencil,1);

maxR = max(stencil(:,2)-stencil(:,1)+1);
basepolyncoeff = zeros(nStencils,maxR,maxR);

r = zeros(nStencils,1);

for s = 1:nStencils
   r(s) = stencil(s,2) - stencil(s,1) + 1;
   basepolyncoeff(s,1:r(s),1:r(s)) = basePolynCoeff(stencil(s,:),centerShift);
end

% when smoothIndShift = 1, cell center type reconstruction
% when smoothIndShift = 2, vertx center type reconstruction
sigma = classicSmoothnessIndSwitch(uBar,stencil,targetCell,basepolyncoeff,smoothIndShift);

eta_bias = zeros(nStencils,1);

hatWgt = linWgt;

eps0 = 1.0;

for s = 1:nStencils
   hatWgt(s) = linWgt(s) / (sigma(s) + eps0*dx^2)^r(s) ...
                         * (eps0*dx / (sigma(s)+eps0*dx))^eta_bias(s);
end

nonlinWgt = hatWgt / sum(hatWgt);

ru = 0;

for s = 1:nStencils
   left  = stencil(s,1) + targetCell;
   right = stencil(s,2) + targetCell;

   uBarStencil = uBar(left:right);

   ru = ru + nonlinWgt(s)*polynEval(hatX,basepolyncoeff,s,uBarStencil,r(s)); 

end

end

% Calculate base polynomial coefficients with a shift from center 
function [sol] = basePolynCoeff(stencil,center)
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

function [ru] = polynEval(p,coeff,s,uBarStencil,r)

   % p is relative coordinate (x-x0)/dx

   ru = 0.0;

   for k=1:r
       for i = 1:r
           ru = ru + uBarStencil(k) * coeff(s,k,i)*p.^(i-1);
       end
   end

end

%{
 {function [val] = polynEval(hatX, basePolynCoeff, uBarStencil, stencil)
 {    % hatX is relative coordinate x/dx
 {
 {    r = stencil(2) - stencil(1) + 1;
 {
 {    val = zeros(length(hatX));
 {    for k=1:r
 {        for p = 1:r
 {            val = val + uBarStencil(k) * basePolynCoeff(k,p)*hatX.^(p-1);
 {        end
 {    end
 {end
 {
 %}
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

function [sigma] = classicSmoothnessIndSwitch(uBar, stencil, ic, basepolyncoeff, s)

    % Create Jiang-Shu Smoothness indicator
    % Shift integral domain with intgShift
    % intgShift = 1 for right side cell
    % inteShift = -1 for left side cell

    switch s
        case 1
            sigma = classicSmoothnessInd(uBar,stencil,ic,basepolyncoeff,1);
        case 2
            sigma = classicSmoothnessInd(uBar,stencil,ic,basepolyncoeff,1);
            sigma = sigma + classicSmoothnessInd(uBar,stencil,ic,basepolyncoeff,-1);
        otherwise
            disp("Invalid type for smoothness indicator"); 
    end

end

function [sigma] = classicSmoothnessInd(uBar, stencil, ic, basepolyncoeff,intgShift)
    gaussPt = [ -sqrt(5 + 2*sqrt(10/7))/3 -sqrt(5 - 2*sqrt(10/7))/3 ...
                0 sqrt(5 - 2*sqrt(10/7))/3 sqrt(5 + 2*sqrt(10/7))/3 ];
    gaussWt = [ (322 - 13*sqrt(70))/900 (322 + 13*sqrt(70))/900 ...
                128/225 (322 + 13*sqrt(70))/900 (322 - 13*sqrt(70))/900 ];
    
    gaussPt = (gaussPt + intgShift)/2;
    %gaussPt = (gaussPt - 1)/2

    gaussWt = gaussWt/2;

    nStencil = size(stencil,1);
    sigma = zeros(nStencil,1);
    
    for s=1:nStencil
        r = stencil(s,2) - stencil(s,1) + 1;
        uBarStencil = uBar(stencil(s,1)+ic:stencil(s,2)+ic);
        %coeff = basePolynCoeff(stencil(s,:),center);
        coeff = reshape(basepolyncoeff(s,1:r,1:r),r,r);

        for ell=1:r-1 %min(2,r-1)
            val = polynEvalDer(ell, gaussPt, coeff, uBarStencil, stencil(s,:));
            sigma(s) = sigma(s) + (val.*val)*gaussWt';
        end    
    end
end
