function f=fact(x);
%> The file <fact.m> calculates the factorials fact(x) of 
%> a matrix with integer entries, sets infinity for negative values.
%> Call: fact(X)
%> Input: X = vector or matrix, integer values
%> Output: a vector or matrix of same size.
%> Compare the standard gamma function of MATLAB. 
%> We want to avoid all NaN outputs!!!
%>
%> © Goran Lindblad - gli@theophys.kth.se

maxx=max(max(x));	y=zeros(size(x));

% Iterate the faculty for elements of values 1,2,...,maxx.

	for n=0:maxx-1

	z=~(x - maxx + n);	y=y+z; 	y=(maxx-n)*y;

	end 

z=~x; % Now deal with the zero arguments, answer is 1.

y=y+z;

% Now deal with all the remaining matrix elements  -
% the arguments are negative, the answer should be infinity -
% set each matrix element separately

w=~y; [w1,w2]=find(w); nn=length(w1);

	for m=1:nn

	y(w1(m),w2(m))=inf;

	end

f=y;


