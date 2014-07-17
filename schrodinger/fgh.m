function varargout = fgh(F,I,M,N)
%> <fgh.m> finds the eigenvalues and eigenfunctions of a 1D Schrodinger 
%> Hamiltonian with a potential defined by the string F. 
%> It uses the Fourier Grid Hamiltonian method, see the reference
%> Marston & Balint-Kurti, J Chem Phys 91(6); 3571 (1989)
%> In this realization we have NOT used the FFT!
%> Call: E = fgh(F,I,M,N)
%> Input: F = a string 'potential(X)' defining the potential
%> 		I = [x1, x2] the interval on the X - axis
%> 		M = an even number of lattice points.
%>		N = the number of eigenvalues to keep
%> Output: A column vector of sorted eigenvalues.
%> Call: [E,V] = fgh(F,I,M,N)
%> Output: eigevalues E and eigenvektors V as N row vectors, each
%> 		of length M+1. 

% Make sure M is even
n2=round(M/2); n0 = 2*n2;  n1=n0+1;
x1=I(1); x2=I(2); 
dx=(x2-x1)/n0;
X=linspace(x1,x2,n1);

%%% evaluate the potential, make a diagonal matrix

pot=eval(F);

%% Setting a large value at the end points to set the eigenstates
%% to zero there
bbb=1.e+3 + max(abs(pot)); 
rr=length(pot); pot(1)=bbb/dx; pot(rr)=bbb/dx;
%%%%%%%

potl=diag(pot);

%%% find the kinetic part of Hamiltonian 

nn1=1:n1; nn2=1:n2;
tt=(nn2*pi/n1/dx).^2 ;
xx1=2*pi*(nn1-1)/n1;
kin=4*tt*cos(nn2'*xx1)/n1;
kin=toeplitz(kin);

%%% add them up 

ham=potl+kin; 


if nargout <= 1

ee=eig(ham);
ee=sort(ee);
varargout={ee(1:N)};


elseif nargout == 2

[vv,ee]= eig(ham);
e1=diag(ee);
[e2,i2]=sort(e1);
e3=[e2(1:N)];
i3=i2(1:N);
v3=vv(:,i3)';

varargout(1) = {e3};
varargout(2) = {v3};

end
