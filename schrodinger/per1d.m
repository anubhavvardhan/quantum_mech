
%> The file <per1d.m> provides a menu for choosing between various 
%> problems for 1D periodic potentials.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear; close;

disp('> Welcome to <per1d>!');
axis('off')
axis([0 1 0 1])
t1=title('Periodic potentials in 1D');
%set(t1,'FontName','courier');
set(t1,'FontSize',18);
hold on

ww = {' '
		' These programs find the Bloch solutions of the Schroedinger equation      ' 
		' for some choices of periodic potentials in 1D and displays the      ' 
		' allowed and forbidden bands.' 
		' ' 
		' The methods used are: ' 
		' (1) Numerical integration using the Numerov algorithm, finding the' 
		' transmission amplitude and then solving for the Bloch waves.' 
		' (2) Solving for the energy eigenvalues for each Bloch momentum' 
		' using matrix approximation.' 
		' (3) Calculating the transfer matrix for discrete models. '};
	

text0([.15 .45  .75   .4],ww); 

rbutt([0.15  0.29  0.35  .06],'Transmission','close;transm')
rbutt([0.15  0.22  0.35  .06],'Bloch solutions','close;kp')
rbutt([0.55  0.29  0.35  .06],'Discrete periodic model','close;trix')
bbutt([.55   .22   .17   .06],'MAIN MENU','close,start')
bbutt([.73  .22   .17   .06],'QUIT','close');


