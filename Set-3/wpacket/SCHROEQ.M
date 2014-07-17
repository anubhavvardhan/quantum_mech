%	Syntax for SCHROEQ.DLL and SCHROEQ.MEXW32 usage:
%			Psi2 = SCHROEQ(Psi1, V, Par);
%
%	INPUT ARGUMENTS:
%	Psi1 = wave function at time t as an N-dim. complex vector; for
%		instance A*exp(-x^2/(4*a^2))*exp(i*k*x); x=-20:0.01:20;
%		a=1; k = C*sqrt(E*Mrel), where E is the mean energy
%		of the particle, Mrel is the mass of the particle relative 
%       to the electron mass m_e=5.11e5 eV, C = sqrt(2*m_e)/(h_bar*c)= 
%       0.5123 in units eV^(1/2)/A; is the electron mass, h_bar= 
%       0.65821 eV*fsec the Planck's constant and c=2998 A/femtosec 
%       (fsec) the light velocity. 
%
%	 V   =	potential energy in eV as a N-dim. real vector (any shape).
%	
%	 Par = [dt,  xWidth,  Mrel,  N], i.e. a real 4-vector, where: 
%		dt = differential time interval in femtoseconds [typically 0.0001];
%		xWidth = width of X-axis interval in Angstrom [typically 40];
%		Mrel = particle-mass/electron-mass ratio [typically 1];
%		N = number of times the incremental solution is computed [e.g. 50]. 
%
%	 OUTPUT ARGUMENT:
%	 Psi2 = the wave function at time t+N*dt [fsec].
%
%	----------------------------------------------------------------------- 
%	This help file describes the usage of the DLL routine SCHROEQ.DLL,
%   called by Matlab versions 5.x and 6.x, or of the MEX routine 
%   SCHROEQ.MEXW32 called by Matlab versions > 7.1. 
%	These routines solve by incremental iterations the 1-dim. 
%   Schroedinger equation: 
%
%	i*h_bar*d Psi(x,t)/dt = -(h_bar^2/2*m)d^2 Psi(x,t)/dx^2 + V(x)*Psi(x,t),	
%
%	where h_bar is the Planck's-constant/2*pi, for a wave amplitude Psi(x,t)
%	moving on an interval I of the x-axis under the action of a potential V(x).
%	The wave undergoes reflection at the ends of I.
% 
%	Actually, the differential equation solved by the routine is:
%
%	i*d Psi(x,t)/dt = -(C1/Mrel)*d^2Psi(x,t)/dx^2 + C2*V(x)*Psi(x,t)
%
%	with  C1 = 5.7884 A^2/fsec,  C2 =1.5193 1/(eV*fsec). This is the
%	form taken by the equation when: the time unit is 1 femtosecond (fsec),
%	the space unit is 1 Angstrom (A), the potential unit is 1 electron-Volt
%	(eV) and Mrel is the particle-mass/electron-mass ratio. Since dt is small
%	but not infinitesimal, the solution is computed using the Caley's form
%	for the unitary elementary translation U(dt): Psi(t+dt)=U(dt)*Psi(t),
%	i.e., with U(dt)=[1-i*H*(dt/2)]/[1+i*H*(dt/2)], where H = Hamiltonian/h_bar. 	
%
%	Author: R.Nobili - Dept. of Physics of Padova University, November 1999.
%   Last revision November 2010.
%