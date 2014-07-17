function f=jacobi(a,b,n);
%> <jacobi.m> calculates the coefficients of the powers of 
%> (1-X)/2 in the Jacobi polynomial of order n. 
%> Defined eg in Edmonds, Angular Momentum in Quantum Mechanics or
%> HMF Chapter 22.
%> Call: jacobi(a,b,n) where a,b real > - 1, n non-neg integer.
%> Output: a matrix of size (n+1)x(n+1); the m:th row consists of the 
%> coefficients of the polynomical of order m-1;
%> The order of the coefficients is from n to 0.
%> For evaluation write something like 
%> z=jacobi(a,b,n),  
%> y=polyval(z(n,:),0.5*(1-X)); X argument which can be a matrix
%> giving a matrix of size = size(X);
%> or y = evalpol(z,0.5*(1-X)) when X is scalar or a row vector
%> giving a matrix of size (n+1)*length(X);

% GL 961123

f=[]; %  Start on function matrix
w=[]; %  start on function vector
a1=1;
a2=1;
nn=[0:n]';
g=gamma(nn+a+1)./fact(nn)/gamma(a+1);% ones(size(nn));

	aa=-nn; % Parameters in hypergeometric series = aa,bb,cc
	bb=a+b+1+nn;
	cc=a+1;

	if n==0
	f=[g,f]; 

	else

	f=[g,f];
		
		for r=1:n
		
		a1=a1*(a+r);
		a2=a2*r;
		g=g.*(aa+r-1).*(bb+r-1)/(cc+r-1)/r;
		f=[f,g];
			
		end



	end

f=fliplr(f);

%%% © Göran Lindblad 1996
