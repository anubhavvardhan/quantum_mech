function f = trf2(X,V,E);
%> This file calculates the transfer matrix for 1D Schroedinger DE
%> for an interval with a linear potential 
%> Call: trf2(x,V,E) 
%> Input: X=[x1, x2] is the interval as a row 2-vector, x1 < x2
%> V = potential values at endpoints =[V(x1), V(x2)];
%> E = energy lattice, row vector
%> Output: a column 4-vector [u;Du;v;Dv] of values at x2
%> where (u,v) are the fundamental solutions defined by
%> u(x1) = Dv(x1) = 1, Du(x1) = v(x1) = 0.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 971222

BB=50; %% Breakpoint, change here

q=0;

if V(2) < V(1)

V=fliplr(V);

q=1;

elseif V(2)==V(1)

f=trf1(X,V(1),E);

return

end

dX=X(2)-X(1); dV=V(2)-V(1);

nn=(2*dV/dX)^(-1/3);

V1= 2*V*nn^2; E1= 2*E*nn^2; 
uniE=ones(size(E));

X10=X/nn;

X1=X10(1) + V1(1) - E1;
X2=X10(2) + V1(1) - E1;

% max(X2), min(X1) 

if max(X2) > BB &  min(X1) <= 3

disp('> The range of values is too large!');

f=[];

else

if max(X2) > BB & min(X1) > 3

y1=aias(X1); y2=aias(X2);

Z1=2*X1.^1.5/3; Z2=2*X2.^1.5/3;

EE=exp(Z1-Z2);

wron1=uniE./pi;


w=[  y2(1,:).*y1(4,:).*EE - y2(3,:).*y1(2,:)./EE;
          y2(2,:).*y1(4,:).*EE - y2(4,:).*y1(2,:)./EE;
        - y2(1,:).*y1(3,:).*EE + y2(3,:).*y1(1,:)./EE;
        - y2(2,:).*y1(3,:).*EE + y2(4,:).*y1(1,:)./EE]./([1;1;1;1]*wron1);

elseif max(X2) <= BB

y1=[airy(0,X1);airy(1,X1);airy(2,X1);airy(3,X1)];
y2=[airy(0,X2);airy(1,X2);airy(2,X2);airy(3,X2)];

wron1=uniE./pi;
 
 w=real([  y2(1,:).*y1(4,:) - y2(3,:).*y1(2,:);
           y2(2,:).*y1(4,:) - y2(4,:).*y1(2,:);
         - y2(1,:).*y1(3,:) + y2(3,:).*y1(1,:);
         - y2(2,:).*y1(3,:) + y2(4,:).*y1(1,:)]./([1;1;1;1]*wron1));
end 
      
if q==0

f=w;

elseif q==1

f=-[w(4,:);w(2,:);w(3,:);w(1,:)];

end

end 

