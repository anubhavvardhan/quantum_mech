function  f=radial(n,x);
%> This file calculates the radial functions of the hydrogen atom 
%> in suitable atomic units, for a range of angular momenta.
%> Call: radial(n,x), 
%> Input: n = positive integer, x = column vector.
%> Output: a matrix of n rows corresponding to  L = 0,...,n-1,
%> columns of size x.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

y=2*x./n; 
z=exp(-0.5*y);

w=[];

for L=0:n-1

k=2*L+1; p=n-L-1;
norm=2*sqrt(prod(1:n-L-1)/prod(1:n+L))/n^2;
w=[w,norm*(y.^L).*z.*laguerre(k,p,y)];

end
	
f=w;


