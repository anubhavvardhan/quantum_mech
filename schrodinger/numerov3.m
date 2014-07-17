function f=numerov3(x0,x1,N,P,W,E,Y0,DY0);
%> The file <numerov3.m> integrates the 1D Schrodinger equation 
%> for an 'arbitrary' potential using the Numerov method.
%> IT IS JUST A MODIFICATION OF NUMEROV1 WITH SMALLER DEMAND ON MEMORY
%> (by avoiding the creation of large matrices).
%> It also allows the potential to be a function of a parameter.
%> It produces the function values and derivatives at the second end 
%> point as a function of energy (to satisfy boundary condition),
%> it does not keep the solution function.
%> Call: numerov3(x0,x1,N,'pot(X,W)',W,E,Y0,DY0);
%> (x0,x1) = end points of integration interval.
%> N = number of steps (N subintervals of length h).
%> P = 'pot(X,W)' if the potential function is stored in <pot.m>.
%> We can let P be a row vector (size E) for each X!
%> W = parameter as a row vector.
%> E = energy lattice as a row vector, same size as W.
%> Y0,DY0 = starting values in x0 , row vectors length(W).
%> Output: f = [Y;DY] = final values of function and derivative in x1,
%> size = [2,length(E)];
%>
%> © Goran Lindblad - gli@theophys.kth.se

%  GL 961231


% Space lattice >> MUST be before eval(P)!

	X=linspace(x0,x1,N)';	uniX=ones(size(X));

% Integration step 

	h=(x1-x0)/(N-1);

% Potential 
% V=mult*eval(P); % rows - energy, columns - space % this is large!

% Sizes

	K1= ones(size(E));	K0= zeros(size(E));	M=length(E);

% Starting conditions
	V1=strrep(P,'X','X(1)');	
	V2=strrep(P,'X','X(2)');
	V1=eval(V1);	
	V2=eval(V2);
	T1=-h^2 *(E-V1)/6;
	T2=-h^2 *(E-V2)/6; 
	F1=(1 - T1).*Y0;
	F2=(1-T2).*((1+6*T1+6*T1.^2).*Y0 + h*(1+2*T1+1.2*T1.^2).*DY0 );

% The Numerov iteration
%%%%%%%%%%%%%%%%%%%%%%%

		for n=2:N-1 

		Vn=strrep(P,'X','X(n)');
		Vn=eval(Vn);
		U = 72./(6+h^2*(E-Vn)) -10;
		F=U.*F2 - F1;
		F1=F2; F2=F;

		end

%%%%%%%%%%%%%%%%%%%%%%%
	VN=strrep(P,'X','X(N)');
	VN=eval(VN);
	VN1=strrep(P,'X','X(N-1)');
	VN1=eval(VN1);
	
	TN=-h^2 *(E-VN)/6;
	TN1=-h^2 *(E-VN1)/6;
% Finding the function and derivative at the point x1:
	XN= F./(K1-TN);
	Y=XN; % this is the function value at right end point.
	XN1= F1./(K1-TN1); % -- in the point before
	DY=((1 + 6*TN + 6*TN.^2).*XN - XN1 )./(1+2*TN+1.2*TN.^2)/h; 
% this is the derivative
%%%% The function numerov3 is the two row vectors Y, DY 
	f=[Y;DY];
	

