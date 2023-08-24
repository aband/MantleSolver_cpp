function [bL,bR] = boundary(a,k,t)

		  %bL = exp(-k*t).*sin(-1-a*t);
        %bR = exp(-k*t).*sin( 1-a*t);

        bL = 1;
		  bR = 0;

