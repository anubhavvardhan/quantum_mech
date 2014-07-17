function L=legpol(n,x);
%> The function file <legpol.m> calculates the Legendre polynomials.  
%> Call: legpol(n,x) 
%> Input: n = non-negative integer, x = row vector.
%> Ouput: matrix dim (n+1) x length(x),
%> the k-th row is the polynomial of order k-1.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125


x0=ones(size(x));
x1=x;
w=[x0;x1];

	if n==0
	L=x0;
	return

	elseif n==1
	L=w;
	return

	end

	for m=2:n

	y=((2*m-1)*x.*x1 - (m-1)*x0)/m;
	x0=x1;
	x1=y;
	w=[w;y];
	
	end

L = w;

%%% © Göran Lindblad 1996
