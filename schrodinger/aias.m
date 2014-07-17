function f=aias(X);
%> This function gives the asymptotics of Airy functions 
%> for positive  arguments, without the exponential factor.
%> Call: aias(X)
%> Input: X a row vector of positive arguments. 
%> For accuracy min(X) > 2
%> Output: a matrix of dimension 4 x dim(X) =
%> [ exp(-zet)*airy(0,X); exp(-zet)*airy(1,X);
%>   exp(+zet)*airy(2,X); exp(+zet)*airy(3.X)]
%> 	where zet = (2/3)*(X)^1.5
%> Reference Handbook of Mathematical Functions ch 10
%>
%> © Goran Lindblad - gli@theophys.kth.se

test=min(X);

if test < 3
disp('> Inadmissible argument.');
f=[]

else

% Coefficients in asymptotic expansion in HMF Section 10.4.59 et sec

kk=[1:5]';
cc=gamma(3*kk+.5)./(54.^kk)./gamma(kk+1)./gamma(kk+.5);
dd= - cc.*(6*kk+1)./(6*kk-1);
cc=[1;cc]; dd=[1;dd]; cc=flipud(cc); dd=flipud(dd);

X1=ones(size(X));
zet=2*(X.^1.5)/3;
pp=[polyval(cc,-1./zet);polyval(dd,-1./zet);polyval(cc,1./zet);polyval(dd,1./zet)];

XX=[X.^(-.25)/2; - X.^.25/2; X.^(-.25);  X.^.25];

f = XX.*pp/sqrt(pi);

end



