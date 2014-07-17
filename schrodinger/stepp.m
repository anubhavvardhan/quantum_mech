
%> The file <stepp.m> calculates the transmission coefficient for a 
%> potential step using the Numerov algorithm.
%> The scattering potential is constant V = V0 > 0 in x < 0, 
%> constant V = 0 in x > 1. 
%> The potential in [0,1] is input as a string like 'cos(pi*X/2)'. 
%> The particle comes in from x = + infinity, we start the integration
%> from the transmitted wave in x = 0, wave number wave2=sqrt(2*(E-V0)), 
%> and integrate to x = 1, wave number wave1=sqrt(2*E), where we adapt the 
%> incoming and reflected waves and renormalize to find the transmission
%> coefficient. This is displayed as a function of E-V. 
%> For comparision the result for a 'pure step' is displayed first.
%> Reference: Messiah Ch III.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961126

clear, close; disp('> Welcome to <stepp>');

txt={' TRANSMISSION COEFFICIENT FOR A 1D POTENTIAL STEP.' 
' ' 
' This file calculates the energy-dependence of the transmission' 
' coefficent for the 1D Schrodinger equation using the Numerov ' 
' algorithm. ' 
' ' 
' The scattering potential V(X) is a potential step modified by a ' 
' function P(X) of your choice in the interval (0,1) which ' 
' interpolates between the values of 1 and 0. ' 
' There is also a choice of an overall multiplier V of the ' 
' potential.' 
' ' 
' The absolute value of the transmission coefficient is displayed ' 
' as a function of the energy variable E - V '};

tt1=text0([.15 .35 .75  .45], txt);
cbutt(.75 , .01), delete(tt1),

%%% First we do the step potential V = 1 for x < 0, V = 0 for x > 1.
% In this case there is no free parameter in the problem. The analytic
% form of the transmission coefficient is known

wave2=linspace(0,ceil(sqrt(6)))+eps;
Y=wave2.^2/2;
wave1=sqrt(2*(Y+1));
TC = 4*wave1.*wave2./(wave1+wave2).^2;

l1=plot(wave2,TC);
set(l1,'LineWidth',2);
xlabel(' The wave number of the transmitted particle ');
title(' Transmission coefficient as a function of wave number ');
t1=text(0.3*wave2(100), 0.75, 'Transmission coefficient of a single step potential ');
set(t1,'FontSize',12);
cbutt(.75,.01);

q1=2;

    while q1 > 0    
    
    if q1==2
    
str='200';
tt2=text0([.15 .45 .75 .05],'Choose the multiplier of the step, V > 0 ');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .45 .15  .05],str);
multiplier=eval(ee);delete(tt2),

str='500';
tt2=text0([.15 .45 .75 .05],'Choose the number of points in space grid ');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .45 .15  .05],str);
N=eval(ee);delete(tt2);

M=200;
energy=2*multiplier;
wave2=linspace(eps ,sqrt(2*energy),M);
xy=[ min(wave2) max(wave2)  0  1.1];

Y=0.5*wave2.^2; % Energy above V
wave1=sqrt(2*(Y+multiplier));
TC = 4*wave1.*wave2./(wave1+wave2).^2;
W=[TC];

	q1=1;
	
	end 

str='1-X';
tt2=text0([.15 .45 .75 .1],'Choose the interpolating function P(X)');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.45 .45 .45  .05],str);
P=ee;delete(tt2);

XX=linspace(0,1); 
XX1=ones(size(XX)); XX0=zeros(size(XX));
X=XX;
pot=eval(P);
l1=plot([XX-1,XX,XX+1],[multiplier*XX1, multiplier*pot, XX0]);
set(l1,'LineWidth',3);
title('This is your potential')
rbutt([.75 .01 .15 .05], 'WAIT', ''), drawnow,

h=1/(N-1); X=[0:h:1]';
pot=eval(P);
uniX=ones(size(X));
V=multiplier*pot;
Y1= ones(size(Y)); Y0= zeros(size(Y));
X1=Y1;	X2=-wave2*i;

%%%%%%% NUMEROV INTEGRATION %%%%%%%%%%%%%
ww=numerov1(0,1,N,P,multiplier,Y+multiplier,X1,X2);

DD = ww(1,:) + i*ww(2,:)./wave1;
NN = ww(1,:) - i*ww(2,:)./wave1;

TT=2./DD;% this is the transmitted amplitude
RR=NN./DD;% this is the reflected amplitude
W=[W;wave2.*abs(TT).^2./wave1];
error=trapz(wave2,Y1- wave2.*abs(TT).^2./wave1 - abs(RR).^2);
disp('> The deviation from unitarity is characterized by the following error:');
disp(error);

cbutt(.75 , .01);
l1=plot(wave2,W);axis(xy);
set(l1,'LineWidth',2);
title(' The transmission coefficients compared ');
xlabel(' Wave number of transmitted particle ');

%%%%%%%%%% INPUT %%%%%%%%%%%
rbutt([.15 .01  .15  .05],'New P(X)','uiresume, obutt, q1=1;'); 
rbutt([.3  .01   .15  .05],'- multiplier','uiresume, obutt, q1=2;'); 
bbutt([.45  .01 .15  .05],'BACK','uiresume, q1=3;');
bbutt([.6  .01   .15  .05],'MAIN MENU','uiresume, q1=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=0;'); 
uiwait;

		if q1==5
	
		q1=0;
	
		elseif q1==4
	
		start;return;	

		elseif q1==3
	
		scatt;return;
	
		end

	end
	
disp('> Type <stepp> to do this again!');


