%> This is <bound1d.m)!
%> This file calculates bound state eigenvalues/functions for a 
%> Schroedinger particle in 1D using the Fourier Grid Hamiltonian
%> method in the function file <fgh.m>. This means: a discrete lattice
%> approximation of the Hamiltonian and the use of the MATLAB <eig> 
%> routine to find eigenvalues and eigenvectors. 
%> The calculation must be done in an interval sufficiently wide such
%> that the exact eigenfunctions are essentially zero outside. 
%> The WKB approximation is used to get a first estimate of the
%> number of solutions. 
%> Uses <fgh.m>, <wkb.m> ... 
%> 
%> © Goran Lindblad - gli@theophys.kth.se
%> 

% GL 990412

% clear; 


close; format long;  disp('> Welcome to <bound1d>');

exx=[]; q1=2;

while q1 > 0
		
		if q1==2

txt={'BOUND STATE EIGENVALUES AND EIGENFUNCTIONS IN 1D' 
' ' 
' This program called <bound1d> finds a number of bound states of a '
' Schrodinger particle in a 1-D potential defined by a user-specified'
' string. ' 
''
' We use the Fourier Grid Hamiltonian (FGH) method (see the help for'
' references) . The problem is discretized and reduced to a matrix '
' eigenvalue problem which is solved using the MATLAB <eig> '
' routine.'
''
' The FGH method can only be used when the exact eigenfuntions go to'
' zero sufficiently quickly outside a finite interval.'
' The algorithm used here searches for negative energy eigenvalues!'
'  '
' You can make the following choices:'
''
' (1) The sech^2 potential (exactly solvable in some cases).'
' (2) The Morse potential (exactly solvable) '
' (3) A square well potential '
' (4) A harmonic potential (exactly solvable) '
' (5) A double-minimum potential '
' (6) A user-defined potential.'
''
' Choose by pressing the corresponding button!'
''
''};

N=120; % change here if you need!! 
N1=N+1;


tt1=text0([.15 .2 .75 .7], txt);

gbutt([.3  .01 .1 .05],'#1','close;q0=1;');
gbutt([.4  .01 .1 .05],'#2','close;q0=2;');
gbutt([.5  .01 .1 .05],'#3','close;q0=3;');
gbutt([.6  .01 .1 .05],'#4','close;q0=4;');
gbutt([.7  .01 .1 .05],'#5','close;q0=5;');
gbutt([.8  .01 .1 .05],'#6','close;q0=6;');
uiwait, 

if q0==1

%% When the multiplier is n*(n+1), then the spectrum is 
%% E_n = -1, -4, -9, --- , - n^2; 

x1=-6; x2=6;
pot='- 56*sech(sqrt(2)*X).^2';
exx=-[49,36,25,16,9,4,1,0]';

elseif q0==2
%% Morse potential, suitably normalized and translated

x1=0; x2=15;
pot ='50*(exp(-2*(X-1.5)) - 2*exp(- (X-1.5)))';
nn=0:9;
exx=.5*(nn'+0.5).*(20 - (nn'+0.5))-50;

elseif q0==3

x1=-5; x2=5;
pot ='- 30*(sign(X+2) + sign(-X+2))';

elseif q0==4
%% Harmonic potential, with a constant subtracted

x1=-5; x2=5;
pot ='3*X.^2 - 30';
nn=0:11;
exx=sqrt(6)*(nn'+.5) -30;

elseif q0==5

x1=-5; x2=5;
pot ='.4*(X-3).^2 .* (X+3).^2  - 36';

end  


if q0==6

x1=-5; x2=5;
str2 =' .05*X.^4 - 12';
gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume');
tt3=text0([.15 .2 .2 .05],'Edit the string >>>  ' );

str3=edit1([.35  .2  .55  .05],str2);

pot=str3;
delete(tt3),

end 

q1 = 1


end 

if q1==1

% delete(tt1,tt2);


X=linspace(x1,x2,N1);
xx=X;

% evaluate potential;
yy=eval(pot); 
miny=min(yy);

if miny >= 0

disp('> V(X) should have some negative part!');

bound1d; 

return

end 

% Estimate of number of bound states
% Stubbe, J Math Phys 31(5), 1177 (1990)
% y0 = .5*(abs(yy)-yy);
% NB=ceil((x2-x1)*sum(sqrt(2*y0))/pi/N1);

% To make a WKB estimate of the eigenvalues, we set Emax > 0 to catch 
% also weakly bounded solutions which get a positive energy from 
% numerical error.
ew=wkb(xx,yy,0.1,300)';
NN=length(ew);

% Plot the potential
l1=plot(xx,yy); hold on,
axis([x1,x2,miny,5]),
set(l1,'LineWidth',2),
title('This is the potential defined by the string'),
xlabel(sprintf('The number of bound states guessed to be %g',NN)),

if NN==0

disp('> V(X) probably has no bound states!');

bound1d; 

return

end

disp(sprintf('> Found the following %g eigenvalues in WKB approximation: ew = ',NN));
disp(ew);

if isempty(exx)==0

edd = diff(mad(exx',ew'))';

disp('> The exact bound state eigenvalues are exx = ');
disp(exx);
disp('> The errors in WKB eigenvalues = ');
disp(edd);

end

cbutt(.75,.01); 

%> Plot eigenvalues on top of potential
plot(xx,ew*ones(size(xx)),'r --'),
title('WKB eigenvalues superimposed on potential');
% xlabel('Press the button to continue!');

% Calculate the FGH eigenvalues and eigenfunctions
[ef,vv]=fgh(pot,[x1,x2],N,NN);

NN=length(ef);
disp(sprintf('> Found the following %g eigenvalues in FGH approximation: ef = ',NN));

disp(ef);

if isempty(exx)==0
edd = diff(mad(exx',ef'))';

disp('> The errors in FGH eigenvalues = ');
disp(edd);

end 


cbutt(.75,.01), 
plot(xx,ef*ones(size(xx)),'k'),
title('FGH and WKB eigenvalues superimposed on potential'),
xlabel('Press the button to continue!'),
cbutt(.75,.01), 
hold off,

plot(xx,vv), axis tight, xy=axis;
title('FGH eigenvectors'),
xlabel('Press a button to continue!');

bbutt([.6   .01  .15  .05],'MORE','uiresume, q1=1;'); 
gbutt([.75  .01  .15  .05],'CONTINUE','uiresume, q1=0;'); 
uiwait,

if q1==1
obutt, 
rbutt([.75  .01  .15  .05],'WAIT...',''); 
drawnow,

for iter=1:NN

plot(xx,vv(iter,:)), axis(xy),
title('FGH eigenvectors'),
pause(1), 


end 

else 

end 


end %% of q1==1

rbutt([.30  .01  .15  .05],'EDIT V(X)','uiresume, obutt, q1=1;');
rbutt([.45  .01  .15  .05],'REPEAT','uiresume, close, q1=2;');
bbutt([.60   .01  .15  .05],'BACK','uiresume, q1=3;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=0;'); 
uiwait;
exx=[];

	if q1==3

	boundst; return;
	
	elseif q1==1

% x1=-5; x2=5;
str2 = pot; 
gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume');
tt3=text0([.15 .2 .2 .05],'Edit the string >>>  ' );

str3=edit1([.35  .2  .55  .05],str2);

pot=str3;
delete(tt3),
	
	end	

	end
	


disp('> Type <bound1d> to do this problem again. Bye!');

%%% © Goran Lindblad 1996-99
