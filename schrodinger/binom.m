function f=binom(x,y)
%> <binom(x,y> calculates the binomial coefficient
%> (x over y) = gamma(x+1)/gamma(y+1)/gamma(x-y+1)
%> x,y  matrices of same dimension, or one scalar, 
%> x-1,y-1,x-y-1 must have strictly positive elements
%> If you use it for integer arguments, round it afterwards.
%> For integers use nchoosek(x,y)
if max(x) < 11

f=gamma(x+1)./gamma(y+1)./gamma(x-y+1);

else

f= exp(gammaln(x+1)-gammaln(y+1)-gammaln(x-y+1));

end
