
%> The file <scatt3d.m> provides a menu for choosing between various
%> scattering problems in 2D and 3D.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear; close;

disp('> Welcome to <scatt3d>!');
axis('off')
axis([0 1 0 1])
t1=title('Scattering problems in 2D and 3D');
%set(t1,'FontName','courier');
set(t1,'FontSize',18);
hold on

ww = {''
' These programs find the partial wave scattering solutions of the ' 
' Schroedinger equation for some central potentials. ' 
'   ' 
' The methods used are: numerical integration using the Numerov ' 
' algorithm, evaluating a few cases where the solutions  to the ' 
' differential equation are known, and the Born approximation.  '    
'  ' 
' In addition, the diffraction from an arbitrary number of slits' 
' is found from the analytical solution.' 
' '};


text0([.15 .45  .75  .4],ww); 


rbutt([.15  .29  .35  .06],'Partial waves in 3D','close;pw3d;')
rbutt([.15  .22  .35  .06],'Hard sphere in 3D','close;hard;')
rbutt([.15  .15  .35  .06],'Square well in 2D','close;pw2d;')
rbutt([.55  .29  .35  .06],'Square well in 3D ','close;well;')
rbutt([.55  .22  .35  .06],'Diffraction from slits','close;diffrac;')
bbutt([.55  .15  .17  .06],'MAIN MENU','close,start')
bbutt([.73  .15  .17  .06],'QUIT','close')


