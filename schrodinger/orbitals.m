
%> The file <orbitals.m> displays hydrogen orbitals indexed by 
%> quantum numbers  [N,L,M]. 
%> It calculates the value (real), with sign, of the 
%> H orbital indexed by integers N,L,M (0 <= M <= L <= N-1)
%> for phi = 0 and draws a color map of the amplitude in the 
%> first quadrant   0 <= r cos theta , r sin theta <= S
%> The length scale S is estimated from N,L in Bohr radius units.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 960226
clear; close; disp('> Welcome to <orbitals>');
str='[1,0,0]';

txt={' HYDROGEN ORBITALS ' 
' ' 
' This file displays the geometric form of hydrogen orbitals.' 
' It draws a color map of the amplitude in the first quadrant,' 
' with a non-linear scale for the amplitude in order to obtain'
' a better contrast for the different orbitals.' 
'' 
' The wave function is indexed by the standard quantum numbers' 
' [N,L,M], N > L >= M .' 
'' 
''
''
' Set the values [N,L,M] in the box at right >>  ' };

tt1=text0([.15 .45  .75 .4],txt);
% title(' Hydrgen orbitals'),
q=1;

ee = str;
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.7 .5  .2  .05],ee);
Q=eval(ee);
delete(tt1);

if isempty(Q)
Q=[1,0,0];
elseif length(Q)~=3
Q=[max(Q),0,0]
elseif Q(1)<=Q(2);
Q=[Q(2)+1,Q(2),0];
elseif Q(2) < abs(Q(3))
Q=[Q(1),Q(2),Q(2)]
end 

		while q == 1

Q(3)=abs(Q(3));
P='The first quadrant amplitude for N = Q(1), L = Q(2), M = Q(3), units of Bohr radius';
P=sprintf(strrep(P,'Q(1)','%g'),Q(1)); 
P=sprintf(strrep(P,'Q(2)','%g'),Q(2)); 
P=sprintf(strrep(P,'Q(3)','%g'),Q(3)); 
% disp(P);


figure(gcf);
N=Q(1);L=Q(2);M=abs(Q(3));
scale =1.1*(2*N*(N+1) - 0.5*L*(L+1) + 0.2*M*(M+1) + 4);
x=linspace(0,scale,50);
y=x;
w=hydrogen(N,L,M,x,y);
mm=max(max(abs(w)));
w=30*w./mm; %scaling can be changed to suit taste.
w=asinh(w); % makes the scale for amplitude logarithmic, with sign!
surf(x,y,w); axis('tight');  axis('equal');
colormap(jet);
view(2);% view([1,-1,3]);
shading interp;%colorbar; % xlabel(P);
% title('Hydrogen orbital amplitude, length scale in units of Bohr radius');
title(P), xlabel('Radial coordinate'), ylabel('Z coordinate');

gbutt([.75 .01 .15 .05], 'PASS', 'uiresume; q=2'),
tt2=text0([.3 .85 .6 .05],'Choose new values for  [N,L,M] >> ');
ee=edit1([.7 .85 .2  .05],ee);
Q=eval(ee);delete(tt2);

if isempty(Q)
Q=[1,0,0];
elseif length(Q)~=3
Q=[max(Q),0,0]
elseif Q(1)<=Q(2);
Q=[Q(2)+1,Q(2),0];
elseif Q(2) < abs(Q(3))
Q=[Q(1),Q(2),Q(2)]
end 

end 

rbutt([.3  .01  .15  .05],'REPEAT ','uiresume; q=1;'); 
bbutt([.45  .01  .15  .05],'BACK','close; q=2;');
bbutt([.6  .01  .15  .05],'MAIN MENU','close;q=3;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q=0;'); 
uiwait;
	
	if q==1 
	
	orbitals;
	
	elseif q==2

	clear; hatom; return;
			
	elseif q==3

	clear; start; return;
	
	end		

disp('> Type <orbitals> to choose another orbital!');


