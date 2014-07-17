%> <start.m> is the starting file for the Schrodinger program.
%> This is a set of files which gives access to some rough and ready MATLAB
%> routines for solving various problems in elementary quantum mechanics.
%> CHOOSE YOUR DEFAULTS in your startup.m file, or right here, to get a
%> window of a suitable size and position, and text of suitable size, eg
%> rect = [70 215 565 573];
%> set(0,'defaultfigureposition',rect);
%> set(0,'defaultaxesfontsize',12);
%> set(0,'defaulttextfontsize',12);
%> You may also have to set the fontsize in <text0.m>, and in the buttons.
%> 
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear; close; global Q1;

% disp('> Welcome to Schrodinger!');
axis('off'); axis([0 1 0 1]);
t1=title('Welcome to Schrodinger!');
set(t1,'FontSize',24); 

ww = {''
' This set of MATLAB programs solves the Schroedinger equation, '
' in 1D, for a choice of potentials and boundary conditions. '
' The methods used are: numerical integration and adaption to boundary  '
' conditions, solving the time-dependent equation by iteration, matrix '
' approximation using a discrete basis, expansion in special functions.'
'               '
' There are a number of nonstandard routines for finding eigenvalues '
' from the numerical integrations. Read the manual and the m-files to '
' learn the basics and the technical tricks. In many of the examples '
' you must know the relevant energy scales in order to set values of '
' the parameters which will give interesting results. '};

xy=[.15 .45 .75 .4];

tt=text0(xy,ww);

rbutt([.15  .38  .35  .06],'Numerical integration','close;integrat')
rbutt([.15  .31  .35  .06],'Particle in a box','close; inbox;')
rbutt([.15  .24  .35  .06],'Wave packets','close;wavepac')
rbutt([.15  .17  .35  .06],'Bound states','close;boundst')
rbutt([.15  .1   .35  .06],'Periodic potentials','close;per1d')

rbutt([.55  .38  .35  .06],'Scattering in 1D','close;scatt')
rbutt([.55  .31  .35  .06],'Scattering in 2D and 3D','close;scatt3d')
rbutt([.55  .24  .35  .06],'Special functions','close;special')
rbutt([.55  .17  .35  .06],'Hydrogen atom','close;hatom')
bbutt([.55  .1   .35  .06],'QUIT','close')


 

