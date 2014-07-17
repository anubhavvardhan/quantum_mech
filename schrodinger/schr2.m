
%> The file <schr2.m> organizes an interactive search for eigenvalues 
%> and eigenfunctions of bound states for the 1D Schrodinger equation. 
%> The potential chosen from  <pot1.m> (it is zero outside [0,1]). 
%> NOTE that it must have some negative values!
%> Method: shooting from 0 and 1 to an interior point using Numerov algorithm.
%> BCs at 0 and 1 adapted to exponential solutions outside.
%> BCs satisfied using <numerov1.m>, solutions calculated using <numerov2.m>. 
%> ERRORS: Positive values of potential create instability problems.
%> <choice.m> allows a choice of potentials in <pot1.m>.
%>
%> © Goran Lindblad - gli@theophys.kth.se


% GL 961123

%help schr2;
%global PP

clear, close, format long, disp('> Welcome to <schr2>');


txt={' INTERACTIVE SEARCH FOR BOUND STATES AND EIGENFUNCTIONS' 
' ' 
' This program organizes an interactive search for bound states and ' 
' eigenfunctions for a particle in a 1D potential whic is zero outside' 
' [0,1]. The potential is chosen from <pot1.m> and must have some' 
' negative values. The BCs: adaption to exponential decay at X = 0,1. ' 
' Method: shooting from 0 and 1 to an interior point, using the' 
' Numerov algorithm, for a lattice of energy values. This is ' 
' followed by an interactive search for the energy values for' 
' which the relevant Wronski determinant is zero. ' 
' ' 
' Now choose your potential!'};

tt1=text0([.15 .5 .75 .35], txt);
cbutt(.75, .01), delete(tt1),

%%% CHOICE OF POTENTIAL 

q2=2;

	while q2 > 0
 
 	if q2==2 	

choice; 

	end
	
subplot(2,2,3),	
XX=linspace(0,1); % grid of 100 points
l1=plot(XX,pot1(XX));
set(l1,'LineWidth',3),
title('This is your potential V(X), normalized.'),

txt={'Choose a multiplier for V(X), with sign, there must be some negative values!!'};
tt1=text0([.15 .55 .75 .1], txt);

gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='-200';
ee=edit1([.75 .55 .15  .05],str);
multiplier=eval(ee); delete(tt1),

l1=plot(XX,multiplier*pot1(XX));
set(l1,'LineWidth',3);
title('This is your potential, with multiplier.')


%%%%%%%%%%%%%%%

W=[]; %starting vector of eigenvalues
YY=[]; %starting matrix of eigenfunctions
% XX=linspace(0,1); % grid of 100 points
% l1=plot(XX,multiplier*pot1(XX)); set(l1,'LineWidth',3);
% title('This is your potential');
% xlabel('Press the button to continue!');

cbutt(.75,.01),

MIN=min(multiplier*pot1(XX)) ;

	if MIN > - eps
txt={'You need a potential with some negative values!!'};
tt1=text0([.15 .55 .75 .1], txt);
cbutt(.75,.01),
schr2; return,

return,

	end
	
disp('> Note that all energies in this file are negative!');
q2=0;
q1=2; %variable number of eigenvalues  q1-1;

while  q1 > 1 %starts procedure for another eigenvalue

if isempty(W)
disp('> Still no eigenvalues...');
else
disp('> The eigenvalues obtained so far are as follows:');
disp(W); 
end

L=[MIN,0];

q3=1;

while q3 < 2

% Define interval of eigenvalue search, and evenly spaced grid

L1=L(1);L2=L(2);
% figure(gcf);
E=linspace(L1,L2,200); % sufficient number of points for interpolation
E0=zeros(size(E));
E1=ones(size(E));
Z=[]; %=empty vector of function values

%%%%%%%%%%%%%%%%%%%% INTEGRATION %%%%%%%%%%%%%%%%%
% Next the DE is integrated for all energy grid points, 
% starting from both endpoints, using BCs and finally setting 
% identifying the logarithmic derivatives at an interior point
% which here is x=0.6 (may have to be changed for some potentials).
% 

K=sqrt(-2*E); Y=E1; DY=E1.*K;
% ww1=numerov1(0,0.6,300,'pot1(X)',multiplier,E,Y,DY);
% ww2=numerov1(1,0.6,200,'pot1(X)',multiplier,E,Y,-DY);
ww1=numerov1(0,0.52,260,'pot1(X)',multiplier,E,Y,DY);
ww1=real(ww1);
ww2=numerov1(1,0.52,240,'pot1(X)',multiplier,E,Y,-DY);
ww2=real(ww2);

ww=(ww1(1,:).*ww2(2,:) - ww1(2,:).* ww2(1,:));%

subplot(1,1,1), plot(E,[ww;E0]),
xlabel('Choose an interval with clearly separated intersections!');

rbutt([.2  .01  .2  .05],'Zoom in','uiresume, q3=1;'); 
rbutt([.4  .01  .2  .05],'Panorama','uiresume, q3=3;');
rbutt([.6  .01  .15  .05],'Zeros','uiresume, q3=2;');
gbutt([.75  .01  .15 .05] ,'CONTINUE','uiresume, q3=4;');
uiwait,

if q3==4

rbutt([.2 .01  .2   .05],'More eigenvalues','uiresume, q3=4;');
rbutt([.4  .01  .2   .05],'New multiplier','close, q3=5;'); 
rbutt([.6  .01  .15   .05],'New V(X)','close, q3=6;'); 
bbutt([.75  .01  .15  .05],'BACK','close, q3=7;'); 
uiwait;

end

if q3 < 5
obutt,
end 

		if q3==7;

		boundstate; return

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

end  %%of q3


% Use findzero to find the zeros (equal logarithmic derivatives)
			if q1==2
			zz=findzero(E,ww)';
			% disp(zz);
			W=[W; zz]; %adds to the vector of eigenvalues
% 			disp('> New eigenvalues:');
			disp(zz);

			end

if q1==1|q1==2
tt1=text0([.15 .5 .75 .35], 'The new eigenvalues are:');
W=sort(W);
			
if isempty(W)

tt2=text0([.15 .5 .3 .3], 'No eigenvalues found!!');

q1=2;

else

tt2=text0([.15 .5 .3 .3], zz);

end

rbutt([.6   .01  .15  .05],'MORE','uiresume;q1=2;'),
gbutt([.75  .01  .15  .05],'ENOUGH','uiresume;q1=1;'), 
uiwait, delete(tt1,tt2),

end %% of q1 | q2 clause

			end %% of q1

% Display eigenvalues

if q2==0

POT=multiplier*pot1(XX');
WW=ones(size(XX'))*W';
WWW=[WW,POT];

plot(XX',WWW),
title('Eigenvalues superimposed on potential'),
xlabel('Wait, then press the button to continue!'),
obutt,
rbutt([.75  .01  .15  .05],'WAIT..',''), drawnow;

E=W'; % eigenvalues
E0=zeros(size(E));
E1=ones(size(E));
Z=[]; %=empty vector of function values

% Next the solutions are calculated, also using the Numerov method:

K=sqrt(-2*E)+eps;
Y=E1;
DY=E1.*K;

% ww1=numerov2(0,0.6,300,'pot1(X)',multiplier,E,Y,DY);
% ww1=ww1./(ones(300,1)*ww1(300,:));
% ww2=numerov2(1,0.6,200,'pot1(X)',multiplier,E,Y,-DY);
% ww2=ww2./(ones(200,1)*ww2(200,:));

ww1=numerov2(0,0.52,260,'pot1(X)',multiplier,E,Y,DY);
ww1=ww1./(ones(260,1)*ww1(260,:));
ww2=numerov2(1,0.52,240,'pot1(X)',multiplier,E,Y,-DY);
ww2=ww2./(ones(240,1)*ww2(240,:));

ww=[ww1;flipud(ww2)];

%%% NORMALIZE!
norm1=ww'*ww/500;
m1=ww(1,:)'*ww(1,:) + ww(500,:)'*ww(500,:);
m2=zeros(length(E));
m2(1,:)=K;
m2=cumsum(m2);
m2=m2+m2';
m1=m1./m2;
norm2=norm1+m1;
norm=sqrt(diag(norm2));
norm=diag(norm);
ww=ww*inv(norm);


[n1,n2]=size(ww);

cbutt(.75,.01);
rbutt([.75  .01  .15  .05],'WAIT..',''),
drawnow,

xx=linspace(0,1,n1);

		for n=1:n2

		plot(xx,ww(:,n)); 		title(sprintf('Eigenfunction # %g',n));

		pause(2);

		end		

plot(xx,ww'); title('All calculated eigenfunctions')

rbutt([.2   .01   .2   .05],'New multiplier','uiresume, obutt, q2=1;'); 
rbutt([.4   .01   .2   .05],'New potential','uiresume, q2=2;');
bbutt([.6   .01   .15  .05],'BACK','uiresume, q2=0;'); 
bbutt([.75  .01   .15  .05],'QUIT','close;q2=-1;'); 
uiwait;

end

end

	if q2==0

	boundst;	return;	

	end

disp('> Type <schr2> to do this problem again!');


