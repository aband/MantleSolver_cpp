function [basepolyncoeff] = basePolyn(stencil,centerShift,grid)

% Take local grid as input instead

nStencils = size(stencil,1);

maxR = max(stencil(:,2)-stencil(:,1)+1);
basepolyncoeff = zeros(nStencils,maxR,maxR);

r = zeros(nStencils,1);

cStart = 1-min(min(stencil));

for s = 1:nStencils
   r(s) = stencil(s,2) - stencil(s,1) + 1;
   basepolyncoeff(s,1:r(s),1:r(s)) = polyn(cStart,stencil(s,:),centerShift,grid);
end

end

function [sol] = polyn(cStart,stencil,centerShift,grid)
    r = stencil(2) - stencil(1) + 1;

    M = zeros(r,r);

    % Extract target cell length
    hCenter = grid(cStart+1)-grid(cStart);
    cCenter = grid(cStart) - hCenter*centerShift;

	 for j = 1:r
		  % Get left and right vertx index
		  vLeftInd  = cStart + stencil(1) + j - 1;
		  vRightInd = cStart + stencil(1) + j;

		  xLeft = grid(vLeftInd);
		  xRight = grid(vRightInd);

		  xLeftShift = (xLeft-cCenter)/hCenter;
		  xRightShift = (xRight-cCenter)/hCenter;

		  for p=1:r
				M(j,p) = xRightShift^p/p - xLeftShift^p/p;
		  end
	 end

    hHat = diff(grid)/hCenter;

	 sol = zeros(r,r);

	 for k = 1:r
		  B = zeros(r,1);
		  B(k) = hHat(k);
		  sol(k,:) = M\B;
	 end

end
