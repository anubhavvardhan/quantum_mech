
%> The file <besjexp.m> gives a series expansion of a function on [0,1]
%> in Bessel functions J_m (calculated by the MATLAB routine besselj(m,x)).
%> The boundary condition is: u=0 at x=1; The function file <besjz.m>
%> is used to find the roots r of the equation J_m(r)=0;
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear; close; format long;

disp('> Welcome to <besjexp>');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FIXED PARAMETERS, YOU CAN CHANGE HERE
%%%%%%
N = 300;		%% Number of points in numerical integration of coeffs
nn = 7;		%% Number of functions in basis displayed.
maxx = 30;	%% Maximal number of terms in expansion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


m0='0'; q1=3;

while q1 > 0

		if q1 == 3
		
txt={' EXPANSION IN BESSEL FUNCTIONS'
''
' This file <besjexp.m> gives a series expansion of functions on [0,1]' 
' in Bessel functions J_m, using the MATLAB function besselj(m,x)).' 
' The boundary conditions aare: u=0 at x=1; The function file ' 
' <besjz.m> is used to find the roots r of the equation J_m(r)=0, '
' '};

t1=text0([.15 .5 .75 .35], txt);
		
txt={' First you have to choose the order m of the Bessels,'
	' non-negative number'};
	
t2=text0([.15 .5 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume'),
m0=edit1([.75 .5 .15 .05],m0); m=eval(m0); delete(t1,t2);

string= ' Calculating the zeros of the Bessel functions up to order %g!';

string=sprintf(string, maxx);
txt={string}; t2=text0([.15 .5 .75 .15], txt);
rbutt([.75 .01 .15 .05], 'WAIT..', ''),
drawnow;		

zeros=besjz(m,maxx);

q1=2; %% 

bbutt([0.5  0.01  0.2  0.05],'Display basis','uiresume; q2=1;'), 
bbutt([0.7  0.01  0.2  0.05],'Expansion','uiresume; q2=2;'),
uiwait; obutt, delete(t2),
		
		while  q2==1

%% Change here! Just for display.
X=linspace(0,1,200);
zz=zeros(1:nn); % zz=besjz(m,nn); 

ww=besselj(m+1,zz); ww=ww.^2/2; arg=zz*X;
basis=besselj(m,arg)./(sqrt(ww)*ones(size(X)));

plot(X,basis); title('The first Bessel functions in the basis');
xlabel('Press the button to continue!'); xy=axis; 
cbutt(.75,.01);
rbutt([0.75  0.01  0.15  0.05],'WAIT..',''), drawnow;

 			for iter=1:nn

plot(X,basis(iter,:)); axis(xy);
title(sprintf('Eigensolution number %g',iter));pause(1);

			end

rbutt([0.6  0.01  0.15  0.05],'REPEAT','uiresume; q2=1;'); 
gbutt([0.75  0.01   0.15   0.05],'CONTINUE','uiresume; q2=2;');
uiwait;

	if  q2==3

	q1=0; close; disp('> Type <besjexp> to do this again!');return;

	end 
	
		end; %% of q2=1
		
		end  %% of q1=3

		if q1==2
		
ff='X.^3.*(1-X)';
		
txt={' A given function will now be expanded in the basis of Bessel functions' 
	' by numerical integration of the expansion coefficients. ' 
	' We have chosen 300 points in the integration interval [0,1].' 
	' Give a real function defined on [0,1] in the editing window' };
		
t1=text0([.15  .5 .75 .25],txt);		
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume'),
ff=edit1([.3  .5 .6 .05],ff),
delete(t1);

X=linspace(0,1,N); fcn=eval(ff); 
fmax=max(fcn); fmin=min(fcn); dd=fmax-fmin;

q1=1; n0='7';
		
		end %% of q1 = 2 

if q1==1
string=' Give the number of terms  < %g you need in the expansion !';
string=sprintf(string, maxx);
t1=text0([.15  .5 .75 .15],string);
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume'),
n0=edit1([.75  .5 .15 .05],n0);
maxn=eval(n0);

if maxn > maxx
maxn=maxx
end

delete(t1);

l1=plot(X,fcn,'-. r'); title('This is the function you have defined');
set(l1,'LineWidth',2);
xy=[0 1 fmin-0.1*dd, fmax + 0.1*dd];
axis(xy); drawnow;hold on;

zz=zeros(1: maxn); ww=besselj(m+1,zz); 
ww=ww.^2/2; %% Norming constant
arg=zz*X; basis=besselj(m,arg); %% Matrix of eigenfunction values 

fc1=(X.*fcn)';fc1(1)=0.5*fc1(1); fc1(N)=0.5*fc1(N);

coeff=basis*fc1/(N-1); %% Trapezoidal rule integration
coef1=coeff./ww; %% Normalization

disp('> The expansion coefficients are:'); disp(coef1);

disp('> The relative error in norm squared is measured by:');

error=1 - coef1'*coeff*(N-1)/(fcn*fc1) ; disp(error);
expand=coef1'*basis; %% Expanded function
plot(X,expand);hold off; axis(xy);
string='The %g - term expansion compared with the function!';
string=sprintf(string,maxn); 
title(string);
xlabel('Press the button to continue!');
cbutt(.75 ,.01);
plot(X,fcn-expand); title('The error term in the expansion');

end

rbutt([.15  .01  .15  .05],'More terms','uiresume; q1=1;'); 
rbutt([.3  .01   .15  .05],'New f(X)','uiresume; q1=2;');
rbutt([.45  .01  .15  .05],'REPEAT','close;q1=3;'); 
bbutt([.6  .01   .15  .05],'BACK','close;q1=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=5;'); 
uiwait;

if q1 < 3

obutt;

end

	if  q1==5

	close; q1=0; 

	elseif q1==4
	
	special; return;
	
	end
	
	end
	
	disp('> Type <besjexp> to do this again');

%%% © Göran Lindblad 1997
