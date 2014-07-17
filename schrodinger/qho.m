
%> The file <qho.m> gives a demo of certain properties of the 
%> 1D QHO eigenstates, and expansions in a basis of such states.
%>
%> © Goran Lindblad - gli@theophys.kth.se

%  GL 961126

clear; close; disp('> Welcome to <qho>!');
 format long; q=3;

while q > 0

		if q==3

txt={' HARMONIC OSCILLATOR EIGENFUNCTIONS AND EXPANSIONS' 
' ' 
' ' 
' This file displays QHO eigenfunctions and gives a series' 
' expansion of an arbitrary function' 
' ' 
' Do you want to see more eigenstates?' 
' Press the appropriate button!'};

tt1=text0([.15 .5 .75 .35], txt);

N=8;
xmax=round(sqrt(N)+3);
x=linspace(-xmax,xmax);
y=ho(N-1,x); % 

%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3),
plot(x,y(1:4,:)),
title('Eigenfunctions n = 0,1,2,3'),
axis('tight'), xy=axis;

%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4),
plot(x,y(4:8,:)),
title('Eigenfunctions n = 4,5,6,7'),
axis(xy),

rbutt([.6  .01  .15  .05],'YES','uiresume;q0=1;'); 
bbutt([.75 .01  .15 .05],'NO','uiresume;q0=0;'); 
uiwait; delete(tt1),obutt,
		
if q0==1 

tt2=text0([.15 .5 .75 .05], 'Choose the number of eigenstates');
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='15';
ee=edit1([.75 .5 .15  .05],str); delete(tt2),
N=eval(ee);

xmax=round(sqrt(2*N)+2);
x=linspace(-xmax,xmax);
% y=hp(N-1);

	
y=ho(N-1,x); % 

subplot(1,1,1),
plot(x,y),
title(sprintf('The %g first QHO eigenfunctions',N)),
xy=axis, rbutt([.75 .01 .15 .05], 'WAIT..','');
xlabel('The eigenstates are localized essentially inside the range of classical motion');
pause(2);

		for n=1:N

		plot(x,y(n,:));
		axis(xy);
		title(sprintf('The QHO eigenfunction #%g ',n));
		pause(1)

		end
xlabel('Press the button to continue!');
cbutt(.75,.01),

end; % of q0


		q=2;

		end 

		if q==2
		
string1={' EXPANSION IN QHO EIGENFUCNTIONS' 
' ' 
' The expansion of an arbitrary function in the eigenfunctions' 
' of the takes place on the whole real line. In order to define ' 
' a problem which can be handled by numerical integration,' 
' define a function which is essentially different from zero' 
' only on the finite interval, say [-10,10].'};

%%%%%%% INPUT
tt1=text0([.15 .5 .75 .35],string1);

gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
string = 'sin(5*X+1)./(1 + .5*(X+.5).^2)';
ee=edit1([.15 .5 .75  .05],string);
string=ee; delete(tt1);

% M=input('> Give the number of integration points >>> ');
M=500;


X=linspace(-10,10,M);X1=ones(size(X));
% close,
y=eval(string); 
subplot(1,1,1), plot(X,y),
title('This is the function to be expanded in oscillator eigenfunctions'); 
xlabel('Press the button to continue!');


		q=1;
		
		end 
		
		if q==1


gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
%%%%%%% INPUT
string = '10';
tt2=text0([.15 .8 .75 .05],'Give the number of terms in expansion');
ee=edit1([.75 .8 .15  .05],string);
N=eval(ee); delete(tt2),

NN=[1:N]; w=ho(N-1,X);
yy=ones(size(NN'))*y;

%% calculate the expansion coefficients
cc=yy.*w; cc=trapz(X',cc'); % Integration by trapezoidal rule
cc1=cc'*X1; 
expand=cc1.*w; 
w=sum(expand);

delta= sqrt(20*sum((w-y).^2)/M);

disp(sprintf('> The RMS error is %g ', delta));
disp('> The expansion coefficients are:');
disp(cc(1:N)');

plot(X,[w;y]), 
title(sprintf('Expansion up to %g terms',N)),
xlabel(sprintf('The RMS error is %g',delta)),
cbutt(.75,.01),

plot(X,w-y),
title(sprintf('This is the error term, RMS value = %g',delta));
	
	end
	
rbutt([.15  .01  .15  .05],'More terms','uiresume, obutt, q=1;'); 
rbutt([.3   .01  .15  .05],'New f(X)','uiresume,obutt, q=2;');
rbutt([.45  .01  .15  .05],'REPEAT','uiresume, obutt, q=3;'); 
bbutt([.6   .01  .15  .05],'BACK','close, q=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','close, q=0;'); 
uiwait;

		if q==4

		special; return;
		
		end

end 

disp('> Type <qho> to do this again!');



