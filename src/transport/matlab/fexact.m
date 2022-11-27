function [u] = fexact(a,k,x,t)

		  u = exp(-k*t).*sin(x-a*t);
