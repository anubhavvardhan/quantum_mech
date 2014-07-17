function f=ylm(L,x);
%> Call f= ylm(L,x) 
%> Input: L = non-negative integer 
%> x = row vector with values in [-1,1].
%> Output: ylm(L,x) = matrix of dimension (2*L+1) x length(x).
%> The program calculates spherical harmonics for
%> m = -L,..,L,  except phi-dependent factor.
%> Ylm = ylm.*exp(i*m*phi); m=[-L:L]'.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125

mm=[0:L]';
y=sqrt((2*L+1).*fact(L-mm)./fact(L+mm)/4/pi);
y=legendre(L,x).*(y*ones(size(x)));
%size(y)
y1=y(2:L+1,:);
y1=(cos([1:L]'*pi)*ones(size(x))).*y1;
y=flipud(y);
y=[y;y1]; 

f=y;


