function  f=hydrogen(N,L,M,x,y);
%> The file <hydrogen.m> calculates the amplitudes of the eigenfunctions 
%> of the hydrogen atom. The amplitude is renormalized by an extra
%> factor = sqrt(sin(theta)*radius^2) (corresponding to the factor
%> sin(theta)*r^2 in the integration element)
%> Call: hydrogen(N,L,M,x,y) 
%> Input: N,L,M integers, 0 ¾ M ¾ L ¾ N-1
%> x,y = row vectors, Bohr radius units, 
%> Output: a matrix of real amplitudes in the plane \phi = 0, 
%> [x,y] = r[cos(\theta) , sin(\theta)]
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 960326

unix=ones(size(x)); uniy=ones(size(y));

rr=(x.^2)'*uniy + unix'*y.^2 + eps;
rr=sqrt(rr);
cos=(x'*uniy)./rr;
rr=2*rr./N;
k=2*L+1;  p=N-L-1;
norm=2*sqrt(fact(N-L-1)/fact(N+L))/N^2;
mm=exp(-0.5*rr);
mm=mm.*(rr.^L);% 
% mm=mm.*sqrt(rr.*(unix'*y)+eps);%this is the extra normalization
mm=norm*mm.*laguerre(k,p,rr);
mm=mm.*yl(L,M,cos);

f=mm;


