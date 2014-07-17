function  f=findzero(x,y);
%> The file <findzero.m> finds 'all' zeros in a given interval.
%> This is done in a rough and ready way by spline interpolation.
%> Call: w=findzero(x,y); finds zeros of y as a function of x.  
%> Input: x, y  = row vectors of same dimension, x representing
%> equispaced arguments of the independent variable
%> Output: w = row vector of zeros.
%> Source of errors: zeros too close cannot be resolved, and are likely
%> to escape detection altogether!
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

d1=diff(y); % 
d2=sign(d1); d3=diff(d2); d4=find(d3);
d4=[0,d4,length(x)-1]+1; d5=length(d4);

z=[];

	for m=1:d5-1

			if d4(m+1) > d4(m)+2 & y(d4(m))*y(d4(m+1)) <= 0

			xxx=x(d4(m):d4(m+1));
			yyy=y(d4(m):d4(m+1));
			z1=spline(yyy',xxx',0);
			z=[z,z1];

			end

	end


f=z;


