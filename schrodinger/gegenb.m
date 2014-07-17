function f=gegenb(a,n);
%> This file calculates the coefficients of the Gegenbauer 
%> polynomial C(a,n), HMF Chapter 22
%> Call gegenb(a,n),  a > 0 , n >= 0 integer
%> Output: (n+1)x(n+1) matrix with integer entries =
%> the coefficients of the n+1 first Gegenbauer polynomials,
%> i.e. the polynomials up to order n inclusive,
%> the k:th row = coefficients of the k:th polynomial
%> To get the polynomials call evalpol(f,x);

%> © Goran Lindblad - gli@theophys.kth.se

if n < 0
disp('> The argument n of gegenb(a,n) must be non-negative!'),
return,
elseif round(n) - n ~= 0
disp('> The argument n of gegenb(a,b) must be an integer!'),
return,

elseif a <= 0
disp('> The argument a of gegenb(a,n) must be positive!');
return

end 


if n==0
w = 1;
else
w = [0  1; 2*a  0];
end


for m=2:n
y=(2*(m -1 + a)*[w(m,:) 0] - (m+2*a-2)*[0  w(m-1,:)])/m;
w=[zeros(m,1) w; y];
end

f=w;
