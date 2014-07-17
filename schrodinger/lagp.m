function f=lagp(k,p)
%> <lagp.m> calculates lagp(k,p), the coefficients of the Laguerre
%> polynomials L(k,0) ... L(k,p)
%> Call P = lagp(k,p);  
%> input p,k non-negative integers
%> output a (p+1)*(p+1) matrix of coefficients
%> To evaluate the value of the polynomial:
%> Call pol = evalpol(P,x);
%> input  x scalar or row vector 
%> output a (p+1)*length(x) matrix evaluating all the polynomials, 
%> OR 	pol=polyval(P(n,:),x) for each row, 
%> input: x a matrix 
%> output: a matrix of same dimension
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 970110

pp=[0:p]'; 

upp=ones(p+1,1);

mm=[0:p];	mm=fliplr(mm);

umm=ones(1,p+1);

D=fact(pp+k)*cos(mm*pi); 

f=D./(fact(pp*umm - upp*mm).*(upp*(fact(mm+k).*fact(mm))));




