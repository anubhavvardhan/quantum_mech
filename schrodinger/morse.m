% <morse.m> calculates the known analytical formula for the eigenvalues of the 
% Morse oscillator, then finds the numerical values using Numerov integration 
% and automatic search for eigenvalues.
% Morse potential: V(x)=k1*(exp(-2*k2*x)-2*exp(-k2*x)); 
% Eigenvalues: E(n)= 0.5*k2*(n+0.5).*[2.*sqrt(2*k1) - k2.*(n+0.5)]-k1; 
% References:
% Morse: Phys Rev vol 34, page 58, 1929.
% Flugge: Practical quantum mechanics, Problem 70.
% Feagin: Quantum methods with Mathematica, Ch 13.
% Normalization: hbar = mass = 1.
%>
%> © Goran Lindblad - gli@theophys.kth.se


%% GL 961125

% clear;
close; disp('> Welcome to <morse>');
%format long;

q=1; str='[30,1]';

%	while q==1
	
txt={' EIGENVALUES OF MORSE POTENTIAL' 
' ' 
' The file <morse.m> displays the known analytic energy eigenvalues ' 
' of the Morse oscillator potential, then calculates the eigenvalues' 
' by numerical integration and compares the results. ' 
' The Numerov method is used for the integration.'
' ' 
' The 2-parameter family of Morse potentials is defined as follows - ' 
' V(x) = k1*[ exp(-2*k2*x) - 2*exp(-k2*x)]' 
' You now have to choose the POSITIVE parameters  [k1,k2]'
''
' Some choices give a bound state of energy 0, for instance'
' [k1,k2] = [50,4]'};

tt1=text0([.15 .45 .75 .4], txt);
cbutt(.75, .01);
delete(tt1);

	while q==1

tt2=text0([.15 .45 .75 .05],'Set the values [k1,k2] > ');
gbutt([.75  .01 .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.7 .45 .2  .05],str);
str=ee; kk=eval(ee);
k1=kk(1); k2=kk(2); delete(tt2),

x=linspace(-0.8/k2,6/k2,30);
X1=ones(size(x));
y= k1*(exp(-2*k2*x) - 2*exp(-k2*x));

l1=plot(x,y,'b');set(l1,'LineWidth',2),
title('This is the shape of the Morse potential with the chosen parameter values.'),
hold on;

N=floor(0.5 + sqrt(2*k1)/k2);

	if N==0

	tt2=text0([.3 .7 .3 .05], 'There are no bound states');
	cbutt(.75 , .01); delete(tt2);
	disp('> There are no bound states!'); hold off;

	else

NN=[0:N-1]';
ex= 0.5*k2*(NN+0.5).*(2.*sqrt(2*k1) - k2.*(NN+0.5))-k1;
y1=ex*X1;

disp(sprintf('> The number of bound states = %g ',N));

	if N > 0
	
	disp('> The analytic solutions of the eigenvalues are: ex =  '); 
	disp(ex);
		
	end
cbutt(.75, .01);
plot(x,y1,'r --');
title('Eigenvalues (exact = - -, numerical = -) superimposed on potential');
% xlabel('Now calculating the numerical values..');
rbutt([.75  .01  0.15  0.05],'WAIT..',';'); 
drawnow;
%% For the numerical solution we rescale the variables

P='exp(-2*X) - 2*exp(-X)';
mult = k1/k2/k2;

x0=-4; x1=9; x2=0; %% Change here if necessary

kk=linspace(0,1,40*N+100);
%% The number of lattice points can be increased for better precision.

E=-k1*kk.^2/k2/k2;
E0=zeros(size(E));
E1=ones(size(E));
N0=100; %% Increase here for better accuracy.
N1=N0*(k2+1);
N1=min(N1,1000);

Y0=E0; DY0=E1;
 
Y1=E1; DY1=-sqrt(-2*E +eps);

w1=numerov1(x0,x2,N1,P,mult,E,Y0,DY0); 
w2=numerov1(x1,x2,5*N0,P,mult,E,Y1,DY1);

ww=w2(2,:) - w1(2,:).* w2(1,:)./w1(1,:);%e

zz=findzero(kk,ww)';
zz= -k1*zz.^2;
en=flipud(zz);
disp('> The eigenvalues calculated with Numerov algorithm are: en =');
disp(en);
N2=length(en);

cbutt(.75, .01);
	y2=en*X1;
	plot(x,y2,'k');
	hold off; 	
	
	
	if N2 == N

disp('> The differences between exact and Numerov eigenvalues are:');
disp(ex-en);

	else
	
disp('> The numerical method did not find the same number of solutions!');

	end

	end

rbutt([.3   .01  .15  .05],'New [k1,k2]','uiresume, obutt, q=1;'); 
bbutt([.45  .01  .15  .05],'BACK','close; q=2;');
bbutt([.6   .01  .15  .05],'MAIN MENU','close; q=3;');
bbutt([.75  .01  .15  .05],'QUIT','close; q=4;');
uiwait;

		end

		if q==2

		boundst;	return;
		
		elseif q==3
		
		start; return;
		
		end
		

	
disp('> Type <morse> to do this again!');
	
