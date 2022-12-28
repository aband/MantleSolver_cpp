% Retrieve coarse mesh from refined mesh from last time step

function [coarseUBar] = coarseMesh(childGrid,childUBar,loc,refineLevel) 

N = length(loc); 

locStart = find(loc,1,'first');
addN     = sum(loc);
locEnd   = locStart+addN; 

coarseUBar = zeros(length(loc),1);

if addN ==0
    coarseUBar = childUBar;
else
    if locStart == 1
        refineStart = 1;
    else
        refineStart = locStart;
        coarseUBar(1:refineStart-1) = childUBar(1:locStart-1);
    end

    if locStart == N 
        refineEnd = N+2;
    else
        refineEnd = refineStart + addN*refineLevel;
        coarseUBar(locEnd:end) = childUBar(refineEnd:end);
    end

    % Recursively average out refined cell values to create coarse cell values
    uBarPreAve = childUBar(refineStart:refineEnd-1);
    uBarPostAve = coarseUBar(locStart:locEnd-1);

    assert(refineLevel*length(uBarPostAve) == length(uBarPreAve));

    for pres = 1:length(uBarPostAve)
        locNow = locStart + pres -1;
        refineSum = 0;
        for posts = (refineLevel*(pres-1)+1):refineLevel*pres
            refineSum = refineSum + uBarPreAve(posts);
        end
        uBarPostAve(pres) = refineSum/refineLevel;
    end

    coarseUBar(locStart:locEnd-1) = uBarPostAve;

end

% End of function
end
