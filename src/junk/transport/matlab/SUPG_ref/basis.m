function [v] = basis(i,p,M,h,type)

v = zeros(size(p));

K = size(p);
KK = K(2);
if type == 1
    for n = 1:KK
    x = p(n);
       if ((x>=(i-2)*h) && (x<=(i-1)*h))
            v(n) = (x-(i-2)*h)/h;
        elseif ((x>=(i-1)*h) && (x<=i*h))
            v(n) = 1 - (x-(i-1)*h)/h;
        end
    end
else if type == 2
    for n = 1:KK
    x = p(n);
       if ((x>=(i-2)*h) && (x<=(i-1)*h))
            v(n) = 1/h;
        elseif ((x>=(i-1)*h) && (x<=i*h))
            v(n) = -1/h;
        end
    end

end

end
