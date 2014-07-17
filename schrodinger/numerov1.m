function f=numerov1(x0,x1,N,P,mult,E,Y0,DY0);
%> The file <numerov1.m> integrates the 1D Schrodinger equation 
%> for an 'arbitrary' potential using the Numerov method.
%> It produces the function values and derivatives at the second end 
%> point as a function of energy (to satisfy boundary condition),
%> it does not keep the solution function.
%> Call: numerov1(x0,x1,N,'pot(X)',M,E,Y0,DY0);
%> Input: (x0,x1) = end points of integration interval.
%> N = number of steps (N subintervals of length h).
%> P = 'pot(X)' if the potential function is stored in <pot.m>.
%> M = multiplying factor for potential, incl sign.
%> E = energy lattice as a row vector.
%> Y0,DY0 = starting values in x0 , row vectors length(E). 
%> Output: f = [Y;DY] = final values of function and derivative in x1,
%> size = [2,length(E)];
%>
%> © Goran Lindblad - gli@theophys.kth.se

%  GL 961125
% last modified 990111


% Integration step 
	h=(x1-x0)/(N-1);
	
		 if h==0
		 
		 f=[Y0;DY0];
		 
		 return
		 
		 end
%% Space lattice >> MUST be before eval(P)!
	X=linspace(x0,x1,N);
%% Potential 
	V=mult*eval(P);
%% Sizes
	uniX=ones(size(X));
	K1= ones(size(E));
	[E1,V1]=meshgrid(E,V);
	
%	M=length(E);
%% The two functions of the Numerov scheme
%	T=-h^2 *(uniX*E-V*K1)/6; 
	
	T=-h^2 *(E1-V1)/6;
	U=(2 + 10*T)./(1 - T);
	
%% Starting conditions
	T1=T(1,:);
%%	T2=T(2,:); T3=T(3,:);
	F1=(1 - T1).*Y0;
 	F2=(1-T(2,:)).*((1+6*T1+6*T1.^2).*Y0 + h*(1+2*T1+1.2*T1.^2).*DY0 );
	
%% Alternative starting values
%% Q Gonzalez & Thompson, Computers in Physics 11(5), 514 (1997)
%% F2=(1-T2).*((1- .5*T3 + 3.5*T1 - 4*T1.*T3).*Y0 + h*(1-T3).*DY0)./(1-3*T2+8*T2.*T3);
%%%%%%%%%%%%%%%%%%%%%%%
%% The Numerov iteration
	for n=2:N-1 

		F=U(n,:).*F2 - F1;
		F1=F2; F2=F;

	end
%%%%%%%%%%%%%%%%%%%%%%%
	clear U;
	TN=T(N,:);
%% Finding the function and derivative at the point x1:
	XN= F./(K1-TN);
	Y=XN; % this is the function value at right end point.
	XN1= F1./(K1-T(N-1,:)); % -- in the point before
%%%% DX=(XN - XN1 + 6*T(N,:).*XN)/h; 

	DY=((1 + 6*TN + 6*TN.^2).*XN - XN1 )./(1+2*TN+1.2*TN.^2)/h; 
%%%% this is the derivative

%%%% The function numerov1 is the two row vectors Y, DY 
	f=[Y;DY];
	
