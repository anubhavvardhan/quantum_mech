function f=yl(L,M,x);
%> Call f=yl(L,M,x)
%> Input: 0 <= M <= L integers, x = matrix, values in [-1,1];
%> Output: yl(L,M,x) = matrix size(x).
%> Calculates the spherical harmonic YLM(cos theta) = yl(L,M,x)
%> except phi-dependent factor (phi = 0).
%> Ylm = yl(L,M,x).*exp(i*M*phi));
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 960326

norm=((-1).^M).*sqrt((2*L+1).*fact(L-M)./(fact(L+M)*4*pi));
f=norm*legf(L,M,x);


