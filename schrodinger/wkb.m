function W = wkb(X,Y,Emax,N);
%> This file finds the semiclassial WKB energy eigenvalues.
%> Call: f = wkb(X,Y,Emax,N)
%> Input: X = lattice in a space interval, Y=pot(X), both row vectors;
%> 	Emax = maximal energy, can be set = 0 or a little bit larger to
%> get a weakly bound state which ends up with a positive energy
%> due to the errors in method and calculation. Emin is automatically 
%> set to min(Y)! 
%>   N = number of points in energy lattice.
%> Output: W = row vector of eigenvalues in increasing order.
%> This procedure should be used with caution unless there is a
%> single minimum of the potential.
%> Reference: For the theory see Messiah Ch 6.
%>
%> © Goran Lindblad - gli@theophys.kth.se

x1=ones(size(X));
Emin=min(Y);
ee=linspace(Emin,Emax,N)';
e1=ones(size(ee));
qq=2*(ee*x1-e1*Y);
hh=0.5*(sign(qq)+1);
qqq=hh.*qq; % this is zero in the forbidden regions.
rqq=sqrt(qqq); % this is the integrand of the WKB formula.
ww=trapz(X,rqq')/pi; % integration by trapezoidal rule.

wmax=max(ww); nmax=floor(wmax-0.5); 

nn=0.5+[0:nmax]; 

W=interp1(ww,ee,nn);

