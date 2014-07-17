% Determining the full width at half maximum
% 
function [FWHM,start] = fwhm(x,f)
if nargin < 2
    f = x;
    x=1:length(f);
end
f1 = abs(f)-0.5*max(abs(f));
ind = find(f1(1:end-1).*f1(2:end) <=0);
FWHM = x(ind(2))-x(ind(1));
start = ind(1);