
%> The file <scatt.m> provides a menu for choosing between various
%> 1D scattering problems.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear, close,

disp('> Welcome to <scatt>!');
axis('off')
axis([0 1 0 1])
t1=title('Scattering problems in 1D');
%set(t1,'FontName','courier');
set(t1,'FontSize',18);
hold on

ww = {' '	
		' These programs find the scattering solutions of the Schroedinger     ' 
		' equation in 1D  and displays the transmission coefficients. ' 
		' ' 
		' The methods used are: ' 
		' (1) Numerical integration using the Numerov algorithm.' 
		' (2) Evaluating a few cases where the analytical solutions     ' 
		' are known.' 
		' (3) Calculating the transfer matrix as a product of contributions' 
		' from intervals where analytic solutions exist. '    
		'                                                                      '};


text0([.15 .45  .75   .4],ww); 


rbutt([.15  .36  .35  .06],'Single well/barrier','close;transm')
rbutt([.15  .29  .35  .06],'Square well/barrier','close; barrier')
rbutt([.15  .22  .35  .06],'Step potential','close;stepp')
rbutt([.15  .15  .35  .06],'Multiple barriers','close;trans1')

rbutt([.55  .36  .35  .06],'Random barriers','close;trans3')

rbutt([.55  .29  .35  .06],'Biased barriers','close;trans2')

rbutt([.55  .22  .35  .06],'Tunneling resonances','close;restun')
bbutt([.55  .15  .17  .06],'MAIN MENU','close,start')
bbutt([.73  .15  .17  .06],'QUIT','close')



