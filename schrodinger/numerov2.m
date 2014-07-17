function f=numerov2(x0,x1,N,P,mult,E,Y0,DY0);
%>  The file <numerov2.m> integrates the 1D Schrodinger equation 
%>  for an 'arbitrary' potential using the Numerov method.
%>  It produces the function values in the integration interval (x0,x1)
%>  for a row vector of energy values.
%>  Call: numerov2(x0,x1,N,'pot(X)',M,E,Y0,DY0);
%>  Input: (x0,x1) = integration interval.
%>  N = number of steps (N subintervals of length h).
%>  P = 'pot(X)' if the potential function is stored in <pot.m>.
%>  M = multiplying factor for potential, incl sign.
%>  E = energy lattice as a row vector.
%>  X0,X1 = starting values in N = x(1),x(2), row vectors size(E).
%>  Output: matrix of dimension (N, length(E));
%>
%> © Goran Lindblad - gli@theophys.kth.se

%  GL 961125

% Integration step 
	h=(x1-x0)/(N-1);
% Space lattice >> MUST be before eval(P)!
	X=linspace(x0,x1,N)';
% Potential 
	V1=mult*eval(P);
% Sizes
	uniX=ones(size(X));
	E1= ones(size(E));	
	M=length(E);

% The two functions of the Numerov scheme
	T=-h^2 *(uniX*E-V1*E1)/6; 
	U=(2 + 10*T)./(1 - T);
	T1=T(1,:);
% Starting conditions
	F11=(1 - T1).*Y0;
	F12=(1 - T(2,:)).*((1+6*T1+6*T1.^2).*Y0 + h*(1+2*T1+1.2*T1.^2).*DY0 );

% Start the matrix of function values in (x0,x1)
	Y1=[F11;F12]; 

% The Numerov iteration
%%%%%%%%%%%%%%%%%%%%%%%
	for n=2:N-1 

		F1=U(n,:).*F12 - F11;
		F11=F12; F12=F1;
		YY= F1./(E1-T(n,:));
		Y1=[Y1;YY];

 	end
%%%%%%%%%%%%%%%%%%%%%%%

f=Y1;

%%% 
