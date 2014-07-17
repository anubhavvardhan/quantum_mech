function [FWHM,nl,nr,n] = findFWHM(x, fx)

f = abs(fx);
[m, n] = max(f);		% max value & index
%FWHM = interp1(fx(end:-1:n), x(end:-1:n), m/2, 'spline') - interp1(fx(1:n), x(1:n), m/2, 'spline');
ind = find(f>=m/2);	        % indicies where I>=max(I)/2
nl = min(ind);			
nr = max(ind);			

%	Linear interpolate x positions
xl = (x(nl)-x(nl-1))*(m/2-f(nl-1))/(f(nl)-f(nl-1)) + x(nl-1);
xr = (x(nr)-x(nr-1))*(m/2-f(nr-1))/(f(nr)-f(nr-1)) + x(nr-1);

FWHM = abs(xr-xl);