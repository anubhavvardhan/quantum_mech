
%> The file <boundst.m> provides a menu of programs calculatating the
%> energies and (in some cases) the eigenstates of some 1D or 3D attractive 
%> potentials. The solution exists as an analytic formula in a few cases,
%> and it can be compared to the numerical values, in other cases there is
%> only a numerical solution. 
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear, close, 
disp('> Welcome to <boundst> !');



axis('off'), axis([0 1 0 1]),
t1=title('Bound states in 1D and 3D');
set(t1,'FontSize',18);
hold on


ww = {' '
		' This set of programs solves for bound state energy levels and  '
		' eigenstates in attractive 1D and 3D potentials. The potential '
		' goes to zero at infinity.'
		' ' 
		' Various methods are used: numerical integration using the Numerov  ' 
		' algorithm and adaption to the boundary conditions, the WKB ' 
		' approximation, numerical search finding the eigenvalues in some   ' 
		' models where the solutions of the DE are known analytically.      ' 
		' In some cases we can compare the exact analytical solutions  ' 
		' with the numerical approximations.  '};

text0([.15 .45 .75 .4] , ww); 

rbutt([.15  .36  .35  .06],'Interactive method in 1D','close;schr2;')
rbutt([.15  .29  .35  .06],'Automated search in 1D','close; schr3;')
rbutt([.15  .22  .35  .06],'Morse potential','close;morse;')
rbutt([.15  .15  .35  .06],'FGH method in 1D','close; bound1d;')
% rbutt([.15  .15  .35  .06],'Chemical bonding, 1D model','close;chembond;')
rbutt([.55  .36  .35  .06],'Spherical square well','close;wellbd;')
rbutt([.55  .29  .35  .06],'3D central potential','close;bound3d')
rbutt([.55  .22  .35  .06],'1D quantum wells','close,well1')
rbutt([.55  .15  .35  .06],'2D quantum wells','close,well2')
bbutt([.55  .08  .17  .06],'MAIN MENU','close;start')
bbutt([.73  .08  .17  .06],'QUIT','close')



%%% © G–ran Lindblad 1996
