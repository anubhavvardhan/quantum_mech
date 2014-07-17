function H = ho(n,x)
%> This file evaluates the QHO eigenfunctions of order up to n.
%> Call: ho(n,x)
%> Input:  n = non-negative integer, x a row vector.
%> Output: a matrix of n+1 columns, rows of dim(x).
%> Units m = omega = hbar = 1. 
%> Uses <hp.m>, <fact.m> and <evalpol.m>
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 990107

nn=[0:n]';

w= 1./sqrt(sqrt(pi)*fact(nn).*2.^nn)*exp(-x.^2/2);

H = w.*evalpol(hp(n),x);

