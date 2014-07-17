function f=evalpol(p,x);
%> evalpol.m does polynomial evaluation for a number of polynomials
%> given by a matrix p, where each row represents the coefficients of a 
%> polynomial, in descending powers.  
%> Call evalpol(p,x); 
%> Input: p a matrix, x a row vector
%> Output = matrix, # rows given by p, # columns = length of x
%> (Compare MATLAB routine polyval) 
%>
%> © Goran Lindblad - gli@theophys.kth.se

pps=size(p); pp1=pps(1); nc=pps(2); nn=[0:nc-1]'; x1=size(x);

f=[];

if sum(x1)==2

	for n=1:pp1

		y = filter(1,[1,-x],p(n,:));
		y = y(nc);
		f=[f;y];
	
	end

elseif x1(1) == 1

xx = (ones(nc,1)*x).^(nn*ones(size(x)));

f=fliplr(p)*xx;

else

disp('> ERROR: the second argument in <evalpol> must be a row vector!');

f=[];

end

