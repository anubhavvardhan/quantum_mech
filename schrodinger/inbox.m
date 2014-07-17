
%> <inbox.m> is the menu file for programs finding eigenvalues and 
%> eigenfunctions for a particle in a 1D box.
%>
%> © Goran Lindblad - gli@theophys.kth.se

close; disp('> Welcome to <inbox>');

axis('off'), axis([0 1 0 1]),
t1=title(' A Schroedinger particle in a 1D or 2D box');
set(t1,'FontSize',18); hold on

ww={''
	' This set of programs solves for the eigenvalues, in some cases also'
	' the eigenvectors of a Schroedinger in a box with infinitely high'
	' walls. The boundary are then that the wave function is zero on '
	' the boundary.' 
	' We solve the 1D problem in the interval [0,1] with a choic of '
	' potentials in the interval.  Two methods are used:            ' 
	'    (1) Solving for a finite matrix in a basis of harmonic functions  ' 
	' satisfying the boundary conditions.                                   ' 
	'    (2) Numerical integration with the Numerov algorithm, numerically ' 
	' solving for the energy values where both boundary conditions are     ' 
	' satisfied. ' 
		' '      
	' The 2D case is solved for a well of quadratic or circular shape'
	' but only for V = 0 inside the boundary'
	                                            
	' There are problems in resolving nearly degenerate eigenvalues and    ' 
	' in the numerical stability in the numerical integrations.            ' 
	' Now choose by pressing a button: '};
	



text0([.15 .4 .75 .45] ,ww); 


rbutt([.15  .29  .35  .06],'1D - matrix approximation','close;matrix1;')
rbutt([.15  .22  .35  .06],'1D - eigenvalues only','close;matrix2;')
rbutt([.15 .15  .35  .06],'1D - Interactive search','close; schr1;')
rbutt([.55  .29  .35  .06],'2D rectangular well', 'well4;')
rbutt([.55  .22  .35  .06],'2D square or circular well','close; well3;')

bbutt([.55  .15  .17  .06],'MAIN MENU','close;start')
bbutt([.73  .15  .17  .06],'QUIT','close')




