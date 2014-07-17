function f=legfun(N,M,x);
%> CALL: f = legfun(N,M,x)
%> Calculates the Legendre functions P(M,L) for L = M,...,N;
%> NORMALIZED to an L2 norm = 1 on the interval [-1,1].
%> Input: N > M integers, x row vector.
%> Output: matrix with N-M+1 rows of length = length(x);
%> This program is just to get this matrix quickly, use the 
%> standard Legendre routines of MATLAB for accuracy.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125


	if M==0

	pmm=ones(size(x));


	else
	
	MM=2*M-1:-2:1;
	ffac= prod(MM);
	pmm=ffac*(sqrt(1-x.^2)).^M;

	end
		
		pmm1=x.*pmm*(2*M+1);
		
		if  N < M
		
		f=[];

		elseif N==M

		f=pmm*sqrt((2*M+1)/2/fact(2*M));

		elseif N==M+1

f= [pmm*sqrt((2*M+1)/2/fact(2*M));pmm1*sqrt((2*M+3)/2/fact(2*M+1))];

		else

		w= [pmm*sqrt((2*M+1)/2/fact(2*M));pmm1*sqrt((2*M+3)/2/fact(2*M+1))];

		for r=M+2:N

		pp=(x.*pmm1*(2*r-1) - (r+M-1)*pmm)/(r-M);
		pmm=pmm1;
		pmm1=pp;
		w=[w;pp*sqrt((2*r+1)*fact(r-M)/2/fact(r+M))];

		end

f=w;

		end


	
	
%%% © Göran Lindblad 1996
