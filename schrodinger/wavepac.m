
%> The file <wavepac.m> provides a menu of files desplaying the 
%> evolution of wavepackets. 
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear; disp('> Welcome to <wavepac>!');

axis('off'), axis([0 1 0 1]);
t1=title('Evolution of wave packets in 1D');;
set(t1,'FontSize',18);
hold on
		
ww = {''
		' This set of programs find the evolution of a wave packet solution' 
		' of the time-dependent Schroedinger equation in 1D for a number' 
		' of different potentials. ' 
		' ' 
		' The solution methods used are: ' 
		' (1) iteration using the Visscher algorithm.' 
		' (2) expansion in a finite subset of eigenfunctions.' 
		'          '};


text0([.15 .45 .75 .4],ww); 


rbutt([.15   .29   .35   .06],'Particle in a box - iteration','close,evolut;')
rbutt([.15   .22   .35   .06], 'Particle in a box - expansion','close, matrix1')
rbutt([.55   .29   .35   .06],'Tunneling - expansion','close,tunnel')
bbutt([.55   .22   .17   .06],'MAIN MENU','close,start')
bbutt([.73   .22   .17   .06],'QUIT','close')


