function  f=fresn(X);
%> Call: fresn(X) = the complex Fresnel function C(X) + i S(X)
%> Input: X is a real vector or matrix      
%> Output: fresn(X) is complex of the same dimension. 
%> For definitions and values of constants, see 
%> Mark A. Heald, Math. Comp. 44(170), 459-461 (1985)
%> and HMF Chap 7. 
%> This is a rational approximation with limited accuracy!
%>
%> © Goran Lindblad - gli@theophys.kth.se

SX=sign(X); X=abs(X);% 
aa=[.0241212,.068324,.2363641,.1945161,1];
bb=[.118247,.356681,.978113,1.875721,2.7570246,2.9355041,2];
cc=[.0433995,.1339259,.3460509,.6460117,.7769507,1];
dd=[.13634704,.4205217,1.0917325,1.9840524,2.7196741,2.5129806,sqrt(2)];
Rx=polyval(cc,X)./polyval(dd,X);
Ax=(polyval(aa,X)./polyval(bb,X) -X.^2)*pi/2;

Cx=.5 - Rx.*sin(Ax);Sx=.5 - Rx.*cos(Ax);

f=SX.*(Cx+i*Sx); 


