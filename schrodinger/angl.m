
%> <angl.m> displays the angular dependence of spherical harmonics for 
%> given input value of the angular momentum L (= non-negative integer).
%> It automatically steps from M = 0 to M = L.  The displayed quantity 
%> is the absolute value of the real part of the function.
%> Uses the function file <ylm.m>
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123


disp('> Welcome to <angl>!');

clear; close;q=0; % help angl;

q=1; 


txt={' ANGULAR DEPENDENCE OF SPHERICAL HARMONICS' 
' ' 
' This file displays the angular dependence of spherical harmonics.' 
' Choose the value of L, it then automatically steps from M = 0 to L.' 
' The displayed quantity is the absolute value of the real part of ' 
' the function YLM(theta, phi).' 
' ' 
' Now set the value of L in the box at right!'};

tt1=text0([.15 .45 .75 .35],txt);
str='3';
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
tt2=edit1([.75 .45 .15  .05],str);
L=eval(tt2);
rbutt([.75 .01  .15 .05], 'WAIT..', ''), drawnow,
delete(tt1);

	while q==1
		
	for M=0:L

%%% Plotting
zz=0.8*[-1 1 -1 1 -1 1]; % defines axis i graphics
%% We choose the lattice for the graphics
theta=pi*linspace(0,1,60); phi=2*pi*linspace(0,1,90); 
sph=ylm(L,cos(theta));
sph=sph(L+M+1,:);
dd=abs(sph' * cos(M*phi));
norm=max(max(dd));
dd=dd/norm;
X=dd.*(sin(theta)'*cos(phi));
Y=dd.*(sin(theta)'*sin(phi));
Z=dd.*(cos(theta)'*ones(size(phi)));

mesh(X,Y,Z),axis('off'),axis(zz),hold on,

l1=line([0,1.1],[0,0],[0,0]);set(l1,'LineWidth',2);

l1=line([0,0],[0,1.1],[0,0]);set(l1,'LineWidth',2);

l1=line([0,0],[0,0],[0,1.1]);set(l1,'LineWidth',2);

l1=text(1.1,0,-0.7,sprintf('L  = %g',L));
set(l1,'FontName','palatino'); set(l1,'FontSize',18),

l1=text(1.1,0,-0.9,sprintf('M = %g',M));
set(l1,'FontName','palatino');set(l1,'FontSize',18),

hold off, 	view([1,1  0.5]),
title('Angular dependence of the spherical harmonics'),	pause(1);

		end % of iteration in M

rbutt([.3   .01  .15  .05],'NEW L ','uiresume; q=1;'),
bbutt([.45  .01  .15  .05],'BACK','close; q=2;'),
bbutt([.6   .01  .15  .05],'MAIN MENU','close;q=3;'),
bbutt([.75  .01  .15  .05],'QUIT','close;q=0;'),
uiwait;

		if q==3
	
		start; return;
		
		elseif q==2
		
		hatom; return;
		
		elseif q==1

obutt,		
tt1=text0([.45 .1 .45 .05],'Set new value in box');		
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
tt2=edit1([.75 .1 .15  .05],str);
L=eval(tt2);
rbutt([.75   .01  .15  .05],'WAIT..',''),
delete(tt1,tt2);		

 		end 

		end % of q

disp('> Type <angl> to do this again!');


