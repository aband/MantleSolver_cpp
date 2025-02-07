% A simple adaptive mesh refinement based on the values of u 
% capturing shock positions
% Cheating by assuming single shock position
function [childGrid, childUBar, loc] = refineMesh(parentGrid, parentUBar, threshold, refineLevel,stencilPolyn32,stencilPolyn43)

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

    N = length(parentUBar);

    stencil32All = zeros(N,size(stencil32,1),size(stencil32,2));

    stencil32All(1,:,:) = stencil32L;
    stencil32All(N,:,:) = stencil32R;

    for i = 2:N-1
        stencil32All(i,:,:) = stencil32;
    end

    % Locate shock location
    loc = (parentUBar<(1-threshold)).*(parentUBar>threshold);
    locStart = find(loc,1,'first');
    addN = sum(loc);
    locEnd   = locStart + addN;

    refineStart = 0;
    refineEnd   = 0;

    if addN == 0
        % No mesh refinement
        childGrid = parentGrid;
        childUBar = parentUBar;
    else
        % Refine base grid with the information of shock location 
        addGrid = linspace(parentGrid(locStart),parentGrid(locStart+addN),addN*refineLevel+1);
        childGrid = [parentGrid(1:locStart),addGrid(2:end-1),parentGrid((locStart+addN):end)];

        % Projection from parent uBar to child uBar
        childUBar = zeros(length(childGrid)-1,1);

        if locStart == 1
            refineStart = 1;
        else
            refineStart = locStart;
            childUBar(1:refineStart-1) = parentUBar(1:locStart-1);
        end

        if locStart == N 
            refineEnd = N+2;
        else
            refineEnd = refineStart + addN*refineLevel;
            childUBar(refineEnd:end) = parentUBar(locEnd:end);
        end

        uBarPreProj  = parentUBar(locStart:locEnd-1);
        uBarPostProj = childUBar(refineStart:refineEnd-1);
        GridPreProj  = parentGrid(locStart:locEnd);
        GridPostProj = childGrid(refineStart:refineEnd);

        assert(refineLevel*length(uBarPreProj) == length(uBarPostProj));

		  % Gauss quadrature points and weights
        gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
        gaussWt = [5/18 8/18 5/18]; %half value!

        % Create Projection
        for pres = 1:length(uBarPreProj)
            vertxLPre = GridPreProj(pres);
            vertxRPre = GridPreProj(pres+1);
            locNow    = locStart + pres -1;
            preCenter = (vertxLPre + vertxRPre)/2;
            for posts = (refineLevel*(pres-1)+1):refineLevel*pres

                vertxLPost = GridPostProj(posts); 
                vertxRPost = GridPostProj(posts+1);

                for g = 1:3
                     gPt = vertxLPost + (vertxRPost-vertxLPost)*(gaussPt(g)+1)/2;
                     gHat = (gPt - preCenter)/(vertxRPre-vertxLPre);
                     uBarPostProj(posts) = uBarPostProj(posts) + gaussWt(g)*...
                     multiLWENORefine(reshape(stencilPolyn32(locNow,:,:,:),nStencils32,maxR32,maxR32),...
                                      parentUBar,reshape(stencil32All(locNow,:,:),3,2), linWgt32,...
                                      locNow,gHat,3,parentGrid);
                end

            end
        end

        childUBar(refineStart:refineEnd-1) = uBarPostProj;

    end

% End of function
end 


