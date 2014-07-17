	function [Psi, x, Ek, lambda] =  GSSPCKT(xo, Dx, Ko, relMass, xWidth, pNum)
%
%	[Psi, x, Ek, lambda] = GSSPCKT(xo, Dx, Ko, relMass, xWdth, pNum)
%
%	Generates a (complex) vector of dimension pNum representing a Gaussian
%	wave-packet centered at xo, with positional incertitude Dx and averaged
%	wave-number Ko. The set of pNum points cover an x-axis interval of width
%	xWidth. The length unit is the Angstrom (A).
%  
% 	INPUT:
%	xo = wave-packet center [def.=0]
%	Dx = standard deviation from xo of the probability distribution [def.=0.9 A]
%	Ko = averaged momentum/h_bar (wave-number) [def.= 0 A^-1]
%	relMass = particle mass/electron mass ratio [def.=1].	
%	xWdth = amplitude of the x-axis interval [def.=40 A]
%	pNum = (odd) number of points covering the x-axis interval [def.=1001].
%
%	OUTPUT:
%	Psi = normalized wave function (pNum-dimensional (complex) vector [column])
%	x = set of x-axis positions (pNum-dimensional real vector [column]).
%	Ek = averaged kinetic energy (electron Volt) 
%	lambda=wave-length (Angstrom).
%			

if nargin < 6,
	pNum=1001;
end

if nargin <5,
	xWidth=40;
end
if nargin <4,
	relMass=1;
end
if nargin < 3,
	Ko=0;
end
if nargin < 2
	Dx=0.9864; % canonical for harmonic oscillator
end
if nargin<1,
	xo=0;
end


N=(pNum-1)/2;
x=xWidth*(-N:N)/pNum;
x=x(:);
x2=(x-xo).^2;
Dx2=Dx*Dx;
Psi = exp(-x2/(4*Dx2));
Norm = sum(Psi.*Psi);
Psi=Psi/sqrt(Norm); % normalized wave function

if Ko==0,
	lambda= 1e16; % wave length in Angstrom
else
	lambda= 2*pi/Ko;
end 

% -----------------------------------
% K = sqrt(2*relMass*m_e*E)/h_bar*c = Co*sqrt(relMass*E),
% where Co = sqrt(2*m_e)/(h_bar*c)= 0.5123 1/(eV^(-1/2)*Angstrom)
% h_bar = 0.6582122 eV*femtosec; c = 2997.92458 Angstrom/femtosec (light speed);
% m_e=5.11e5 eV (elecron mass). E is measured in eV. 	  

h_bar=0.6582122; % Planck's constant in eV*femtosec	

po = Ko*h_bar;   % averaged momentum	
Dp = h_bar/(2*Dx);	% momentum indetermination

c=2997.92458;	   % light velocity in A/femtosec
m_e = 5.11e5;	% electron mass in eV.	

Ek =((po*c)^2 + (Dp*c)^2)/(2*relMass*m_e); % averaged kinetic energy
iK=sqrt(-1)*Ko;	% otherways Psi is real
Psi = Psi.*exp(iK*x);
