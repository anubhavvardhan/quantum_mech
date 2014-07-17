% SCHRODINGER TOOLBOX
% Version of April 1999
%==========================================================================
% <schrodinger> is a toolbox of MATLAB files which which gives the neophyte 
% student of Quantum Mechanics an interactive access to some of the basic 
% methods of solving the Schrodinger equation. Read more in the
% <readme.txt> file and in the  M-files listed below!
%
% The files can be run on any computer with the MATLAB program installed 
% with a graphical interface.  They have been tested with MATLAB 5.2 on a 
% PowerMac 7600 and on a PC/Linux. They should run also on version 4.2.
% The speed evidently depends on the platform, but it can be adapted to
% your computer through the choice of the number of lattice points and 
% basis functions. This choice will inevitably also influence the 
% accuracy of the solutions.
% 
% Start the session by executing the file <start.m> then making choices
% by clicking the buttons and inserting values in editing windows.
%
% You may have to change your window size and the font size in some places!
%
%==========================================================================
% Files in the schrodinger folder, F - says it is a function file.
%==========================================================================
% Contents.m - this file
% Readme.txt - Read this file for an introduction.
% start.m - EXECUTE THIS FILE TO GET STARTED
%==========================================================================
% DEMOS - you can execute these independently of <start.m>
%==========================================================================
% angl.m - plots the angular dependence of spherical harmonics.
% barrier.m - finds the transmission coeff for a 'square' barrier/hole.
% besjexp.m - finds expansions in Bessel functions J_m
% bound1d.m - finds bound states for 1D Schrodinger particle
% bound3d.m - finds bound states for central potential defined in <pot2.m>
% boundst.m - menu file for bound state problems
% choice.m - utility for interactive choice of potential from <pot1.m>
% diffrac.m - calculates the diffraction from an array of slits
% evolut.m - displays the evolution of a wave packet using iteration.
% hard.m - displays scattering cross sections for hard sphere in 3D.
% hatom.m - menu displaying some features of hydrogen eigenstates.
% inbox.m - menu file for particle in a box + potential
% kp.m - finds spectrum of Kronig-Penney like model.
% legexp.m - expansions in Legendre polynomials and functions.
% integrat.m - shows effect of instability on numerical integration
% matrix1.m - solves for eigenvalues/functions using matrix approx.
% matrix2.m - solves for eigenvalues only, faster - uses slider input.
% morse.m - calculates the eigenvalues of Morse potential
% orbitals.m - graphic display of H orbitals.
% per1d.m - menu choice of periodic potential problems
% pw2d - scattering for 2D circular potential well
% pw3d.m - partial wave method in 3Dscattering
% qho.m - displays some properties of QHO eigenstates.
% restun.m - finds narrow tunneling resonances.
% scatt.m - a menu for 1D scattering problems.
% scatt3d.m - menu for 3D scattering problems.
% schr1.m - interactive search for eigenvalues/functions.
% schr2.m - interactive search for bound states in 1D.
% schr3.m - automatic search for bound states in 1D.
% special.m - a menu for some special functions.
% stepp.m - scattering on a step potential plus an additional term.
% trans1.m - transmission for a multiple barrier
% trans2.m - transmission for a biased multiple barrier
% trans3.m - transmission for a multiple barrier with random heights
% transm.m - transmission coefficient and Bloch spectrum in 1D.
% trix.m - calculates Bloch solutions in a discrete model
% tunnel.m - scattering and tunneling of wave packet in 1D.
% wavepac.m - menu file for wave packet evolutions
% well.m - partial wave scattering for potential well.
% well1.m - bound states in various 1D quantum wells.
% well2.m - bound state energies in quadratic and circular 2D wells
% well3.m - eigenvalues of an infinite quadratic or circular 2D well
% well4.m - eigenvalues of infinite rectangular 2D well
% wellbd.m - bound states in a spherical 'square' well.
%==========================================================================
% FUNCTION FILES
%==========================================================================
% aias.m - F - asymptotics for Airy functions for positive arguments.
% besjz.m - F - calculates zeros of Bessel functions J_m
% binom.m - F - calculates binomial coefficients
% crop.m - F - utility which crops values outside a finite range.
% evalpol.m - F  - polynomial evaluation
% fact.m - F - factorial for integer-valued matrix inputs.
% fgh - F - Calculates Schrodinger eigenvalues/states using the
%	Fourier Grid Hamiltonian method
% findzero.m - F - finds zeros by spline interpolation.
% fresn.m - F - calculates the complex Fresnel function
% gegenb.m - F - calculates coefficients of Gegenbauer polys
% glue.m - F - glues jumps of given amplitude in a function
% ho.m - F - evaluates QHO eigenfunctions.
% hp.m - F - calculates coefficients of Hermite polynomials.
% hydrogen.m - F - calculates hydrogen eigenstate amplitudes.
% inv4.m - F - inversion for 2x2 matrices represented as 4-vectors
% jacobi.m - F - calculates coefficients of Jacobi polynomials
% lagp.m - F - coefficients of Laguerre polynomials
% laguerre.m - F - evaluates Laguerre polynomials, matrix arguments
% legf.m - F - calculates Legendre functions for a matrix argument
% legfun.m - F - calculates a set of Legendre functions
% legpol.m - F - calculates the values of a set of Legendre polys
% mad.m - F - adds a row vector to a matrix, adjusting the dimensions
% mult4.m - F - multiplication of 2x2 matrices represented as 4-vectors
% numerov1.m - F - integrates Schrodinger eq, keeps boundary values only
% numerov2.m - F - integrates Schrodinger eq, keeps solution
% numerov3.m - F - variant of numerov1, use when memory gets short
% pot1.m - F - list of potentials used by <choice.m>
% radial.m - F - calculates hydrogen radial eigenfunctions
% radnl.m - F -  another calculation of radial solutions
% sde.m - F - defines a Schrodinger 1D DE for ode45 etc.
% transfer.m - F - calculates transfer matrix using Numerov integration.
% trf1.m - F - transfer matrix for an interval, constant potential
% trf2.m - F - transfer matrix for an interval, linear potential
% wkb.m - F - calculates semiclassical WKB eigenvalues
% yl.m - F - calculates spherical harmonics
% ylm.m - F - - " -, other format only.
% zed.m - F - an auxiliary function used to find zeros of Bessels
% zeroai.m - F - calculates zeros of Airy function
%==========================================================================
% GRAPHICS - INPUTS - TEXT
%==========================================================================
% edit1.m - F - input through editing a string
% cbutt.m - F - a green 'continue' button
% bbutt.m - F - blue button
% gbutt.m - F - green button
% obutt.m - F - grey non-functional button
% rbutt.m - F - red button
% slider1.m - F - input using a slider
% suptitle.m - F - useful graphics (by Drea Thomas)
% text0.m - F - default text box, 12 pt
% wbutt.m - F - a red 'wait' button
%==========================================================================
%%% © Goran Lindblad 1999 - gli@theophys.kth.se
