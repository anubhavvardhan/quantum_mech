function f=zed(x);
%> Call: f=zed(x)
%> Input: x in [1,infty]
%> This is the function z(zwt) defined in Handbook of Mathematical
%> Functions sect.  9.5.22 as the inverse of zet(z) defined in
%> HMF sect. 9.3.38-39.
%> Uses spline interpolation to perform the inverse.
%>
%> © Goran Lindblad - gli@theophys.kth.se

maxx=max(x)+2;

maxy=(1/1.5)*maxx.^1.5;

y=linspace(1,maxy,30);

w=(1.5*(sqrt(y.^2 - 1) - acos(1./y))).^(2/3);

% f=interp1(w,y,x,'spline');

f=spline(w,y,x);


