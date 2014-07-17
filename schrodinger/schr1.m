
%> The file <schr1.m> organizes an interactive search for eigenvalues
%> and eigenfunctions for a Schrodinger particle in a 1D box = [0,1]. 
%> The potential chosen from  <pot1.m>.  BCs at 0 and 1  are u=0;
%> Method: shooting from 0 and 1 to an interior point using Numerov algorithm.
%> BCs satisfied using <numerov1.m>, solutions calculated using <numerov2.m>. 
%> ERRORS: Positive values of potential create instability problems.
% <choice.m> allows a choice of potentials in <pot1.m>.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

clear, close, format long, disp('> Welcome to <schr1>!');

txt={' INTERACTIVE SEARCH FOR EIGENVALUES AND EIGENFUNCTIONS' 
' ' 
' This program organizes an interactive search for eigenvalues and ' 
' eigenfunctions for a particle in a 1D box = [0,1]. ' 
' The potential is chosen from  <pot1.m>.  BCs: u(0) = u(1) =  0.' 
' Method: shooting from 0 and 1 to an interior point, using the' 
' Numerov algorithm, for a lattice of energy values. This is ' 
' followed by an interactive search for the energy values for' 
' which the relevant Wronski determinant is zero. ' 
' ' 
' Now choose your potential!'};

tt1=text0([.15 .5 .75 .35], txt);
cbutt(.75, .01), delete(tt1),

q2=2;

while q2 > 0

%%%%%%%%%%%%% CHOICE OF POTENTIAL 
 
 			if q2==2

			choice; 

			end

subplot(2,2,3),	
XX=linspace(0,1); % grid of 100 points
l1=plot(XX,pot1(XX));
set(l1,'LineWidth',3),
title('This is your potential, normalized.'),

txt={'Choose a multiplier for your potential, with sign!'};
tt1=text0([.15 .55 .75 .1], txt);

gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='-200';
ee=edit1([.75 .55 .15 .05],str);
multiplier=eval(ee); delete(tt1),

l1=plot(XX,multiplier*pot1(XX));
set(l1,'LineWidth',3);
title('This is your potential, with multiplier.')

q2=0;
%%%%%%%%%%%%%%%

W=[]; %starting vector of eigenvalues
YY=[]; %starting matrix of eigenfunctions
XX=linspace(0,1); % grid of 100 points
MIN=min(multiplier*pot1(XX)) ;

txt={'Choose the maximal energy for the eigenvalues!'};
tt1=text0([.15 .55 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');

str='200';
ee=edit1([.75 .55 .15  .05],str);
MAX=eval(ee); delete(tt1),

L0=[MIN,MAX];

q1=2; %variable number of eigenvalues  q1-1;

		while  q1 > 1 %starts procedure for another eigenvalue

if isempty(W)
disp('> Still no eigenvalues...');
else
disp('> The eigenvalues obtained so far are as follows:');
disp(W); 
end

L=L0;
		q3=1;

				while q3 < 2

% Define interval of eigenvalue search, and evenly spaced grid

disp('> Note that the energies here can be of either sign');

L1=L(1);L2=L(2);
E=linspace(L1,L2,200); % sufficient no of points for interpolation
E0=zeros(size(E));
E1=ones(size(E));
Z=[]; %=empty vector of function values

%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%
% Next the DE is integrated for all energy grid points, 
% starting from both endpoints, using BCs and finally setting 
% identifying the logarithmic derivatives at an interior point
% which here is x=0.55 (may have to be changed for some potentials).

				Y=E0; DY=E1;

ww1=numerov1(0,0.55,550,'pot1(X)',multiplier,E,Y,DY);
ww2=numerov1(1,0.55,450,'pot1(X)',multiplier,E,Y,-DY);

% ww1=numerov1(0,0.5,500,'pot1(X)',multiplier,E,Y,DY);
% ww2=numerov1(1,0.5,500,'pot1(X)',multiplier,E,Y,-DY);

ww=(ww1(1,:).*ww2(2,:) - ww1(2,:).* ww2(1,:));%
subplot(1,1,1),
plot(E,[ww;E0]);
xlabel('Choose an interval with clearly separated intersections!');

rbutt([.3   .01  .15  .05],'Zoom','uiresume, q3=1;'); 
rbutt([.45  .01  .15  .05],'Panorama','uiresume, q3=3;');
rbutt([.6   .01  .15  .05],'Zeros','uiresume, q3=2;');
gbutt([.75  .01  .15  .05] ,'CONTINUE','uiresume, q3=4;');
uiwait,

if q3==4

rbutt([.2  .01  .2   .05],'More eigenvalues','uiresume, q3=4;');
rbutt([.4  .01  .2   .05],'New multiplier','close, q3=5;'); 
rbutt([.6  .01  .15  .05],'New V(X)','close, q3=6;'); 
bbutt([.75 .01  .15  .05],'BACK','close, q3=7;');
uiwait;

end

		if q3==7;

		inbox; return

		elseif q3==6

		q1=0;q2=2;
		
		elseif q3==5
		
		q1=0; q2=1; 
					
		elseif q3==4
		
		q1=1; q2=0; 

		elseif q3==3
		
		q1=3; 
	
		elseif q3==2
		
		q1=2;
		
		elseif q3==1
		
		[L,Y]=ginput(2);

		end %% of q3 clause
%
		end  % of q3


% Use findzero to find the zeros (equal logarithmic derivatives)

				if q1==2
		
				zz=findzero(E,ww)'; 
				

% 				disp('> New eigenvalues:');
				disp(zz);
				W=[W; zz]; %adds to the vector of eigenvalues

				end;
				
				if q1==1|q1==2

tt1=text0([.15 .5 .75 .35], 'The new eigenvalues are:');
				
				W=sort(W);
% disp(W);
if isempty(W)
tt2=text0([.15 .5 .3 .3], 'No eigenvalues found!!');

q1=2;

else

tt2=text0([.15 .5 .3 .3], zz);

end

obutt,
rbutt([.6   .01  .15  .05],'MORE','uiresume;q1=2;'),
gbutt([.75  .01  .15  .05],'ENOUGH','uiresume;q1=1;'), 
uiwait, delete(tt1,tt2),

		end % q1|q2
		end % of q1
		
if q2==0;

% Display eigenvalues

disp(W);
POT=multiplier*pot1(XX');
WW=ones(size(XX'))*W';
WWW=[WW,POT];
plot(XX',WWW);
title('Eigenvalues superimposed on potential');
xlabel('Wait .. then press the button to continue!');
obutt, rbutt([.75  .01  .15  .05],'WAIT..',''),
drawnow;

E=W'; % eigenvalues
E0=zeros(size(E));
E1=ones(size(E));
Z=[]; %=empty vector of function values

%Next the solutions are calculated, also using Numerov method:

K=sqrt(2*E);

				Y=E0; DY=E1;

ww1=numerov2(0,0.55,550,'pot1(X)',multiplier,E,Y,DY);
ww1=ww1./(ones(550,1)*ww1(550,:));
ww2=numerov2(1,0.55,450,'pot1(X)',multiplier,E,Y,-DY);
ww2=ww2./(ones(450,1)*ww2(450,:));

% ww1=numerov2(0,0.5,500,'pot1(X)',multiplier,E,Y,DY);
% ww1=ww1./(ones(500,1)*ww1(500,:));
% ww2=numerov2(1,0.5,500,'pot1(X)',multiplier,E,Y,-DY);
% ww2=ww2./(ones(500,1)*ww2(500,:));

ww=[ww1;flipud(ww2)];

%%% NORMALIZE!
norm1=ww'*ww/1000;
norm=sqrt(diag(norm1));
norm=diag(norm);
ww=ww*inv(norm);

[n1,n2]=size(ww);

cbutt(.75 , .01);
rbutt([.75  .01  .15  .05],'WAIT..',''),
drawnow;

xx=linspace(0,1,n1);

				for n=1:n2

				plot(xx,ww(:,n));
				title(sprintf('Eigenfunction # %g',n));
				pause(2);

				end

xx=linspace(0,1,n1);
plot(xx,ww');
title('All calculated eigenfunctions');

rbutt([.2   .01   .2   .05],'New multiplier','uiresume, obutt, q2=1;'); 
rbutt([.4   .01   .2   .05],'New potential','close; q2=2;');
bbutt([.6   .01   .15  .05],'BACK','close;q2=0;'); 
bbutt([.75  .01   .15  .05],'QUIT','close;q2=-1;'); 
uiwait;

end; % of q0 clause

end % of q2

	if q2==0

	inbox;	return;	

	end

disp('> Type <schr1> to do this problem again. Bye!');

