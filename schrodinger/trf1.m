function f=trf1(X,V,E);
%> This file calculates the transfer matrix for 1D Schroedinger DE
%> for an interval with a constant potential 
%> Call:  trf1(X,V,E)
%> Input: X=[X1, X2] is the interval as a row 2-vector, X1 < X2
%> V = potential, scalar
%> E = energy lattice, row vector, INCREASING ARGUMENT
%> Output: a column 4-vector [u;Du;v;Dv] at X2
%> where the fundamental solutions (u,v) are defined by
%> u(X1) = Dv(X1) = 1, Du(X1) = v(X1) = 0.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 971222
  
dx=X(2)-X(1); ww=2*(E-V);

pos=(sign(ww+eps)+1)*.5; w1=find(1-pos);w2=find(pos);

k1=sqrt(-ww(w1)+eps); k2=sqrt(ww(w2)+eps);


f=[ cosh(k1*dx)     ,    cos(k2*dx)     ;
    k1.*sinh(k1*dx) ,  - k2.*sin(k2*dx) ;
    sinh(k1*dx)./k1 ,    sin(k2*dx)./k2 ;
    cosh(k1*dx)     ,    cos(k2*dx)     ];
	

