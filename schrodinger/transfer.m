function f=transfer(x0,x1,N,P,mult,E);
%> The file <transfer.m> integrates the 1D Schrodinger equation 
%> for an 'arbitrary' potential using the Numerov method.
%> It produces the transfer matrix from x0 to x1
%> as a function of energy (to satisfy boundary condition),
%> Call: transfer(x0,x1,N,'pot(X)',M,E);
%> Input: (x0,x1) = end points. Interchanging x0,x1 gives inverse matrix.
%> N = number of steps (N subintervals of length h).
%> P = 'pot(X)' if the potential function is stored in <pot.m>.
%> M = multiplying factor for potential, incl sign.
%> E = energy lattice as a row vector.
%> Output: f = [u(x1); Du(x1); v(x1); Dv(x1)] = transfer matrix 
%> Fundamental solutions (u,v) defined by 
%> u(x0)=1, Du(x0)=0, v(x0)=0, Dv(x0)=1;
%> size of f = [4,length(E)]; 
%>
%> © Goran Lindblad - gli@theophys.kth.se

 
%  GL 961125

% Integration step 
h=(x1-x0)/(N-1);
% Space lattice >> MUST be before eval(P)!
X=linspace(x0,x1,N)';
% Potential 
V=mult*eval(P);

% Sizes
uniX=ones(size(X));
E1= ones(size(E));	
E0= zeros(size(E));
M=length(E);

% The two functions of the Numerov scheme
% T=-h^2 *(uniX*E-V*E1)/6; 
% U=(2 + 10*T)./(1 - T);

T1=-h^2 *(E-V(1))/6; %T(1,:);
T2=-h^2 *(E-V(2))/6; %T(1,:);
F1=(1 - T1);
F2=(1-T2).*(1+6*T1+6*T1.^2) ;

G1= E0;
G2=h*(1 - T2).*(1+2*T1+1.2*T1.^2).*E1;

% The Numerov iteration
%%%%%%%%%%%%%%%%%%%%%%%
	for n=2:N-1 
	%U=(12-10*h^2*(E-V(n)))./(6+h^2*(E-V(n)));
	U = 72./(6+h^2*(E-V(n))) -10;
	F=U.*F2 - F1; G=U.*G2 - G1;
	%F=U(n,:).*F2 - F1; G=U(n,:).*G2 - G1;
	F1=F2; F2=F; G1=G2; G2=G;
	end
%%%%%%%%%%%%%%%%%%%%%%%
clear U 

% Finding the function and derivative at the point x1:
TN=-h^2 *(E-V(N))/6; %T(N,:);
TN1=-h^2 *(E-V(N-1))/6; %T(N-1,:);
XN= F./(E1-TN);u=XN; % this is the u function value at right end point.
XN1= F1./(E1-TN1); % -- in the point before
Du=( - XN1 + (1+6*TN+ 6*TN.^2).*XN)./(1+2*TN+1.2*TN.^2)/h; % this is the derivative 

XN= G./(E1-TN);v=XN; % this is the v function value at right end point.
XN1= G1./(E1-TN1); % -- in the point before
Dv=( - XN1 + (1+6*TN+ 6*TN.^2).*XN)./(1+2*TN+1.2*TN.^2)/h; % this is the derivative

%% The function transfer gives  the 4  row vectors u,Du,v,Dv;
%% The transfer matrix is 2x2: [u, v; Du, Dv].%% [u, Du; v, Dv]
%% The determinant f(1,:).*f(4,:) - f(2,:).*f(3,:) = 1 is a test of stability.
%% The inverse is [Dv, -Du; - v, u] => [f(4,:); - f(2,:); - f(3,:);f(1,:)];
%% The product of two matrices is, of course:
%% [f(1)g(1)+f(2)g(3); f(1)g(2)+f(2)g(4); f(3)g(1)+f(4)g(3); f(3)g(2)+f(4)g(4)];

f=[u;Du;v;Dv];


