function [err] = errorLnorm(order,uBar,a,k,x,t,index)

% Three point gauss quadrature rule
gaussPt = [-sqrt(3/5) 0 sqrt(3/5)];
gaussWt = [5/18 8/18 5/18]; %half value!

if order == 0
    % Return L infinity norm
    err = 0;

    vertxL = x(1:end-1); 
    vertxR = x(2:end);
    for g = 1:3
        gPt = vertxL + (vertxR-vertxL)*(gaussPt(g)+1)/2;
        err = err + gaussWt(g)*abs(fexact(a,k,gPt,t)-uBar);
    end

    err = max(err(index));

else
    % Return L-order norm
    err = 0;

    vertxL = x(1:end-1); 
    vertxR = x(2:end);
    for g = 1:3
        gPt = vertxL + (vertxR-vertxL)*(gaussPt(g)+1)/2;
        err = err + gaussWt(g)*abs(fexact(a,k,gPt,t)-uBar).^order;
    end

    err = sum(err(index)).^(1/order);

end

end
