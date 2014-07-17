function F=pot1(X);
%> The function file <pot1.m> defines a number of potentials on [0,1].
%> Using the file <choice.m> a value of the parameter PP is picked. 
%>
%> © Goran Lindblad - gli@theophys.kth.se

%  GL 961125


global PP pot1str

if PP==1
% F=exp(-0.005.*ones(size(X))./((X+eps).*(1-X+eps)));
% 
F= 2*sqrt(X.*(1-X));

return
% else 
% end
elseif PP==2
F = 3.53*( sqrt(X.*(-X+0.4)+ abs(X.*(-X+0.4)))...
 + sqrt((X-0.6).*(-X+1)+ abs((X-0.6).*(-X+1))));

return
 

elseif PP==3
F = 7.05*( sqrt(X.*(-X+0.20)+ abs(X.*(-X+0.20)))...
+ sqrt((X-0.4).*(-X+0.6)+ abs((X-0.4).*(-X+0.6)))...
+ sqrt((X-0.8).*(-X+1)+ abs((X-0.8).*(-X +1))));
% F= exp(-25*(X-0.5).^2);

return
 

elseif PP==4
% F= sin(2*pi*X).^2;
%%%%
F = 2.81*( sqrt(X.*(-X+0.5)+ abs(X.*(-X+0.5)))...
 + sqrt((X-0.65).*(-X+1)+ abs((X-0.65).*(-X+1))));


return
 

elseif PP==5
% F=  sin(3*pi*X).^2
%
% F = 4.7*( sqrt(X.*(-X+0.27)+ abs(X.*(-X+0.27)))...
% + sqrt((X-0.4).*(-X+0.6)+ abs((X-0.4).*(-X+0.6)))...
% + sqrt((X-0.7).*(-X+1)+ abs((X-0.7).*(-X +1))));
%%%%%%%%
F= sech(8*(X-0.5)).^2; 
%  With a factor - 32*n*(n+1) this is exactly solvable!');
%  The eigenvalues are then - 32*k^2 , k= 0,1,..,n-1.'
%  The potential is also reflectionless.

return 

elseif PP==6
% F= sin(3*pi*X);
F = 9.4*( sqrt((X-0.05).*(-X+0.2)+ abs((X-0.05).*(-X+0.2)))...
 + sqrt((X-0.3).*(-X+0.45)+ abs((X-0.3).*(-X+0.45)))...
 + sqrt((X-0.55).*(-X+0.7)+ abs((X-0.55).*(-X+0.7)))...
 + sqrt((X-0.8).*(-X+0.95)+ abs((X-0.8).*(-X+0.95))));

return 

elseif PP==7
% F=(floor(1+X-eps)).*(1-floor(X)).*4.*(X-0.5).^2;
F=4.*(X-0.5).^2;

return 

elseif PP==8
F=110*( X.^2.*(1-X).^2.*(0.5-X));

return

elseif PP==9
%F=(1-floor(X)).*ceil(X).*(0.5*(sign(X-0.65) - sign(X-0.35))+1);
F=ceil(X-0.01).*ceil(0.99-X).*(0.5*(sign(X-0.65) - sign(X-0.35))+1);

return

elseif PP==10
F = 7.05*( sqrt(X.*(-X+0.15)+ abs(X.*(-X+0.15)))...
 + sqrt((X-0.25).*(-X+0.45)+ abs((X-0.25).*(-X+0.45)))...
 + sqrt((X-0.55).*(-X+0.75)+ abs((X-0.55).*(-X+0.75)))...
 + sqrt((X-0.85).*(-X+1)+ abs((X-0.85).*(-X+1))));
% F= 5*sqrt((X+eps).*(1-X+eps)).* sin(4*pi*X).^6;

return

elseif PP==11
F= X;

return

elseif PP==12

F=eval(pot1str);

end

return



