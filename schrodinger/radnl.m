function  f=radnl(N,L,X)
%> This file calculates the radial functions of the hydrogen atom 
%> in suitable atomic units, for a given angular momentum.
%> Call: radnl(N,L,X), 
%> Input: N = positive integer, 0 <= L integer < N, 
%> X = matrix of positive elements
%> Output: a matrix of same size as X.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 9901

if L >= N

disp('> You must choose  0 <= L < N !!');

f=[];


elseif N-round(N) + L-round(L) ~= 0

disp('> You must choose N,L integers')

f=[];


else


norm=2*sqrt(prod(1:N-L-1)/prod(1:N+L))/N^2;
	
f=norm*(2*X/N).^L .*exp(-X/N).*laguerre(2*L+1,N-L-1,2*X/N);

end


