% A simple adaptive mesh refinement based on the values of u 
% capturing shock positions
function [childGrid, childUBar] = meshRefinement(parentGrid, parentUBar, threshold, refineLevel)

    % Locate shock location
    loc = (uBarCurrent<(1-threshold)).*(uBarCurrent>threshold);
    locStart = find(loc,1,'first');
    addN = sum(loc);

    if addN == 0
				childGrid = parentGrid;
	 else
        % Refine base grid with the information of shock location 
        childGrid = [x(1:locStart),linspace(x(locStart),x(locStart+addM),addM*refineLevel+1),x((locStart+addM):end)];



end 
