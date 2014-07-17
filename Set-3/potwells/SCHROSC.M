%
%	THIS HELP FILE DESCRIBES THE USAGE OF THE *.DLL ROUTINE SCHROSC.DLL 
%   FOR MATLAB VERSION 6.5 OR THE MEXFILE SCHROSC.MEXW32. 
%
%	SCHROSC(Psi_1, t_1, Vo, DV, Par) solves the Schroedinger's Eq.
%	for a 1-Dim wave function Psi defined on a X-axis interval of
%	given width, passed by Par(4), whose potential changes periodically
%	according to the law 
%	                   V(t) = Vo + sin(omega*t)*DV. 
%	The computation is carried out from time t_1 to t_2=t_1+cycleNum*dt,
%	as a sequence of 'cycleNum' unitary transformations, each one of which
%	represents the evolution of the system in an elementary time interval
%	dt, which is passed by Par(1); t_1 is the initial time; omega, passed
%	by  Par(2), is the radian frequency of the potential oscillation in
%	radian/femtosec, to be compared with Bohr's radian frequencies
%	omega_ij=(E_i-E_j)/h_bar which are characteristic of energy level
%	transitions; cycleNum, passed by Par(3), is the number of time the
%	computation is carried out at each call.
%	
%	USAGE:   [Psi_2, t_2, V] = SCHROSC(Psi_1, t_1, Vo, DV, Par);
%
%	with:
%	INPUT QUANTITIES:
%	  Psi_1 = input wave function (N-dim, real or complex) at time t_1
%	  t_1 = starting time of computation in femtosec [1e-15 sec]
%	  Vo = unperturbed potential energy (real N-dim)
%	  DV = potential perturbation	(real N-dim) [eV]
%	  Par = [dt, omega, cycleNum, xWidth, m_r]; 5 parameters:
%	  dt = time interval of the elementary evolution	
%	  omega = radian frequency of the perturbation in radian/fsec;  
% 	  cycleNum = number of cycles the elementary evolution is computed.
% 	  xWidth = width of the wave-propagation interval on X-axis in Angstrom [1e-10 m]
%	  m_r = mass relative to electron mass
%
%	OUTPUT QUANTITIES:	 
%	 Psi_2 = the returned wave function computed at time t_2
%	 t_2 = t_1 + cycleNum*dt, the value of time returned after computation;
%	 V = Vo + sin(omega*t_2)*DV, the potential profile after computation. 
%
%	MORE INFORMATION:
%	The routine SCHROSC.DLL, solves by incremental iterations the 1-dim.
%	Schroedinger equation: 
%
%	i*h_bar*d Psi(x,t)/dt = -(h_bar^2/2*m)d^2 Psi(x,t)/dx^2 + V(x,t)*Psi(x,t),	
%
%	where h_bar is the Planck's-constant/2*pi, for a wave amplitude Psi(x,t)
%	moving on an interval I of x-axis under the action of a potential V(x,t).
%	generally varying in time. The wave undergoes reflection at both ends of I.
% 
%	Actually, the differential equation solved by the routine is:
%
%	i*d Psi(x,t)/dt = -(C1/Mrel)*d^2Psi(x,t)/dx^2 + C2*V(x, t)*Psi(x,t)
%
%	with  C1 = 5.7884 A^2/fsec,  C2 =1.5193 1/(eV*fsec). This is the
%	form taken by the equation when: the time unit is 1 femtosecond (fsec),
%	the space unit is 1 Angstrom (A), the potential unit is 1 electron-Volt
%	(eV) and Mrel is the particle-mass/electron-mass ratio. As  dt is small
%	but not infinitesimal, the solution is computed using the Caley's form
%	for the unitary elementary translation U(t, dt): Psi(t+dt)=U(t,dt)*Psi(t),
%	with U(t, dt)=[1-i*H(t)*(dt/2)]/[1+i*H(t)*(dt/2)], where H(t) =
%	Hamiltonian/h_bar at time t. 	
%	The radian frequency omega appearing in the potential interpolation
%	formula can be expressed in eVolt through the relationship
%	E_Plank = h_bar*omega. This may be of help to relate the radian
%	frequency of the potential perturbation with energy levels. 
%		
%	Author: R.Nobili - Dept. of Physics of Padova University, November 2010. 
%