function f=hp(n);
%> Call: hp(n)
%> Input: n non-negative integer 
%> Output: (n+1)x(n+1) matrix with integer entries =
%> the coefficients of the n+1 first Hermite polynomials,
%> i.e. the polynomials up to order n inclusive,
%> the k:th row = coefficients of the k:th polynomial
%> in descending order. 
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

if n < 0
disp('> The argument of <hp> must be non-negative!'),
return,
elseif round(n) - n ~= 0
disp('> The argument of <hp> must be an integer'),
return,
end 

x0=[1];  x1=[2 0]; % starting values.

	if n==0

	w=x0;

	else

	w=[0 x0;x1];

	end


	for m=2:n
	
	y=2*[x1 0]-2*(m-1)*[0 0 x0];
	x0=x1; 	x1=y;
	w=[zeros(m,1) w; x1];

	end 


f=w;

