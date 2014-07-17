function f=laguerre(k,n,x);
%> The function file <laguerre.m> calculates the values of 
%> the Laguerre polynomials.
%> Call: laguerre(k,n,x), 
%> Input: k = real, n = non-negative integer, x = matrix of values.
%> Output: a matrix of same dimension as x.
%> Reference: HMF Ch 22.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961210
x1=ones(size(x));

	if n==0
	w3=x1;

	elseif  n==1

	w3=k+1-x;

	else

w1=x1; w2=k+1-x;

	for r=2:n

	w3=2*w2-w1+((k-1-x).*w2 - (k-1)*w1)/r;

	w1=w2; w2=w3;

	end

end

f=w3;


