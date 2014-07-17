
%> The file <schr3.m) calculates bound state eigenvalues/functions 
%> for a Schroedinger particle in 1D. The potential is nonzero only on
%> a finite symmetric interval [-A,A].  
%> BC: continuity with bound state exponentials outside [-A,A] and
%> continuity of logarithmic derivative at an interior point  X = B.
%> Method: Shooting from -A and A to X = B using the Numerov algoritm
%> for a lattice of energy values. The zeros of the Wronski determinant
%> at  X = B gives the eigenvalues. The search for solutions is automatic.
%> The method cannot resolve nearly degenerate eigenvalues! 
%> It is likely to miss very weakly bound states!
%> Uses files <wkb.m>, <numerov1.m>, <findzero.m> ......
%>
%> © Goran Lindblad - gli@theophys.kth.se

%%%%%%%%%%%% GL 990425

clear; close; format long;  disp('> Welcome to <schr3>');

q1=2;

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 N=400; %% You can change here, # steps in Numerov algorithm
 A=5; %%%  sets interval [-A,A]; other choices below
 M=500; %%% sets number of energy lattice points
 N0=50; %%% defines an offset point of continuity for Wronskian
 potstr =' X.^4/9 - 9'; %%% defines editable string for potential
 
while q1 > 0
		
		if q1==2

txt={' FINDING EIGENVALUES AND EIGENFUNCTIONS IN 1D' 
' ' 
' This program searches for bound states for a Schrodinger particle  ' 
' in a 1D potential which is approximately zero outside a finite '
' symmetric interval [- A,A]. ' 
' It first finds the WKB approximation to the eigenvalues, then '
' better approximations by numerical integration of the Schrodinger'
' equation. '
''
' Boundary conditions: continuity with bound state exponentials at '
' at end points X = - A, A, continuity of logarithmic derivative '
' at an interior point  X = B. '
' Method: Shooting from - A and A to B using the Numerov algorithm.'
' This is done for a lattice of energy values; spline interpolation'
' finds the zeros of the Wronski determinant at X = B, and hence '
' the energy eigenvalues.' 
' You can make the following choices:'
''
' (1) The sech potential (exactly solvable in some cases).'
' (2) A square well potential '
' (3) A harmonic well potential'
' (4) An asymmetric well potential '
' (5) A user-defined potential.'
' .' 
''};

tt1=text0([.15 .25 .75 .65], txt);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gbutt([.3  .01 .06 .06],'#1','close;q0=1;');
gbutt([.4  .01 .06 .06],'#1','close;q0=1;');
gbutt([.5  .01 .06 .06],'#2','close;q0=2;');
gbutt([.6  .01 .06 .06],'#3','close;q0=3;');
gbutt([.7  .01 .06 .06],'#4','close;q0=4;');
gbutt([.8  .01 .06 .06],'#5','close;q0=5;');
uiwait,  %delete(tt1),

if q0==1

A=6;
%% When the multiplier is n*(n+1), then the spectrum is 
%% E_n = -1, -4, -9, --- , - n^2; 

pot='- 42*sech(sqrt(2)*X).^2';
exx=-[36,25,16,9,4,1,0]';

elseif q0==2

A=3;
pot ='- 6*(sign(X+2) + sign(-X+2))';

elseif q0==3
%% Harmonic potential, with a constant subtracted
A=5;
pot ='.2*X.^2 - 5';
nn=0:11;
exx=sqrt(6)*(nn'+.5) -30;

elseif q0==4

A=3;

pot ='- 6*(sign(X+2) + sign(-X+2)).*(1 + .2*X)';

elseif q0==5

A=3;

gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume');
tt3=text0([.15 .2 .2 .05],'Edit the string >>>  ' );

pot=edit1([.35  .2  .55  .05],potstr);
potstr=pot 

delete(tt3),

end  


x=linspace(-A,A,300);
X=x; y=eval(pot);
MIN=min(y);

% Estimate of number of bound states
NB=ceil(2*A*sum(sqrt(2*abs(y)))/pi/300);
% Stubbe, J Math Phys 31(5), 1177 (1990)
 
% subplot(2,2,3),
% l1=plot(x,y);
% axis([-A,A,MIN,1]),
% set(l1,'LineWidth',2),
% title('This is the potential defined in the string <pot>'),
% xlabel(sprintf('The number of bound states guessed to be  %g',NB)),

%% We make a WKB estimate of the eigenvalues,setting Emax > 0 to catch 
% also weakly bounded solutions which get a positive energy from numerical error.
eigen=wkb(x,y,0.1,M)';
NN=length(eigen);
disp(sprintf('> Found %g eigenvalues in WKB approximation = ',NN));
disp(eigen);
% cbutt(.75,.01), % delete(tt1),


l1= plot(x,y);axis([-A,A,MIN,1]),
set(l1,'LineWidth',2), hold on 
plot(x,eigen*ones(size(x)),'r --'),
% axis([-A,A,MIN,1])
title('WKB eigenvalues superimposed on potential');
xlabel('Press the button to continue!');
cbutt(.75,.01), 
clear eigen  

%% Now for the numerical integration, 
% space grid fixed, energy grid and interval variable


% txt={'Choose the size of the energy grid '};
% tt1=text0([.15 .55 .75 .05], txt);
% gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
% str='400';
% ee=edit1([.75 .55 .15  .05],str);
% M=eval(ee); delete(tt1),
		
% 		disp(sprintf('> M = %g',M));
		
		end %% of q1==2
		
		eigen=[];

%% INPUT
str='[xx,0]';
str=sprintf(strrep(str,'xx','%g'),MIN); 

txt={'Choose an energy interval [E0,E1] for eigenvalues, E0, E1 < 0! '};
tt1=text0([.15 .55 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
ee=edit1([.65 .55 .25  .05],str);
Range=eval(ee); delete(tt1),
rbutt([.75  .01  .15  .05],'WAIT..',''),
drawnow;

disp('> Integrating using the Numerov algorithm!');

B=A*N0/N;%%% offset point of continuity for Wronskian

%% The lattice is chosen in the wave number space
waveno=linspace(sqrt(-2*Range(2)),sqrt(-2*Range(1)),M);
K=-0.5*waveno.^2;	K1= ones(size(K));	K0= zeros(size(K));

Y=K1; DY=K1.*waveno; % Starting values in x = +/- 3.

ww1=numerov1(-A,B,N+N0,pot,1,K,Y,DY);
ww2=numerov1(A,B,N-N0,pot,1,K,Y,-DY);

ww=(ww1(1,:).*ww2(2,:) - ww1(2,:).* ww2(1,:));% Wronski determinant

%% Find the roots, using <findzero.m>, this algorithm will 
%% miss states which are marginally bound, having energy 0: 

zz=fliplr(findzero(waveno,ww))';
eigen=-zz.^2/2;

if isempty(eigen)
disp('> Found no eigenvalues!');

else

NN=length(eigen);
disp(sprintf('> Found %g eigenvalues = ',NN)); 
disp(eigen);

l1= plot(x,y);axis([-A,A,MIN,1]),
set(l1,'LineWidth',2),hold on,
l1=plot(x,eigen*ones(size(x)),'k');
% set(l1,'LineWidth',2);
title('Numerical eigenvalues superimposed on potential'),
xlabel('A weakly bound state may be missing!'), 
hold off,
rbutt([.75  .01  .15  .05],'WAIT..',''),drawnow;

% We now calculate the eigenfunctions on the given interval
% using the same integration method, we just restrict the 
% <waveno> variable to the values corresponding to the 
% eigenvalues <eigen> already calculated above, and save the
% function value at each step

% disp('> Busy calculating the eigenfunctions...');

K=eigen';
unie=ones(size(eigen));
waveno=sqrt(-2*K); 
K0= zeros(size(K)); K1= ones(size(K));

Y=K1; DY=K1.*waveno; % Starting values in x = +/- 3.

ww1=numerov2(-A,B,N+N0,pot,1,K,Y,DY);
ww1=ww1./(ones(N+N0,1)*ww1(N+N0,:));
ww2=numerov2(A,B,N-N0,pot,1,K,Y,-DY);
ww2=ww2./(ones(N-N0,1)*ww2(N-N0,:));

YY=[ww1;flipud(ww2)];
%% size(ww1), size(ww2)

clear ww1 ww2

h=2*A/(2*N-1);	XXX=[-A:h:A]; %new x coordinates, 2*N-1 points
% Now we want to normalize the eigenfunctions. Start with diagonal
% elements in order to have same order of magnitude 

norm=[];

		for m=1:NN
		norm(m)=sqrt(h*YY(:,m)'*YY(:,m));%???
		end

size(norm); norm=diag(norm);
YY=YY/norm'; %???

% Now we calculate the normalization matrix including the 
% contribution of the exponential parts outside [-A,A].

D1=h*YY'*YY;
D2=YY(1,:)*YY(1,:)' + YY(2*N,:)*YY(2*N,:)';%???
D3=unie*waveno + waveno'*unie';
D2=D2./D3;	D1=D1+D2;	
% size(D1);	contour(D1) % This should be a unit matrix ideally
norm=diag(sqrt(diag(D1))); YY=YY/norm;
cbutt(.75,.01);

l0=plot(XXX,[YY,zeros(size(XXX'))]);
% axis('off');
title(sprintf('All the %g eigenfunctions',NN));
hold off,

bbutt([.6   .01  .15  .05],'MORE','uiresume, q1=1;'); 
gbutt([.75  .01  .15  .05],'CONTINUE','uiresume, q1=0;'); 
uiwait,

if q1==1

rbutt([.75  .01  .15  .05],'WAIT..',''),
drawnow;

		for n=1:NN;

		plot(XXX,YY(:,n));
		title(sprintf('This is eigenfunction # %g',n));
		xlabel('Check the number of nodes!');
		pause(2)

		end

end

end

% disp('> That is all there is, folks!');

rbutt([.3   .01  .15  .05],'SET [E0,E1]','uiresume, obutt, q1=1;'); 
rbutt([.45  .01  .15  .05],'REPEAT','uiresume, close, q1=2;');
bbutt([.6   .01  .15  .05],'BACK','uiresume, q1=3;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=0;'); 
uiwait;

	if q1==3

	boundst; return;
	
	end	

	end

disp('> Type <schr3> to do this problem again. Bye!');

%%% © Goran Lindblad 1996-99
