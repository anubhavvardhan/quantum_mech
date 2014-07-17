
%> The file <tunnel.m> displays the scattering/tunneling in a square 
%> potential well/barrier combined with a step.
%> We use the analytic form of the solution as given e.g. in Messiah. 
%> The continuous Fourier analysis in terms of scattering states is 
%> discretized, hence replaced by a finite Fourier series. 
%> The energy eigenvalues are known, so we get the unitary evolution.
%> This methods contrasts with the solution by iteration of the evolution.
%> The advantage here is the fast display, the disadvantage is the storing
%> of a matrix of a very high dimension (space lattice x wave number lattice). 
%> WARNING: You will need a good deal of memory to get the graphics good!
%> The amplitudes involves division by wave numbers, when these are near zero
%> things get a bit shaky.
%> The wave packet comes in from the right. 
%> You can set 3 parameters, two energy parameters for the potential,
%> and the initial velocity.
%> A1 = potential in [-5,0], it is zero in [1,5],
%> A2 = height/depth of barrier in [0,1].
%> k0 = velocity of wave packet.
%> Parameters which can be changed in file: 
%> N, M give the number of lattice points, 
%> D defines the width of the wave packet.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear;close; disp('> Welcome to <tunnel>!');

txt={' TUNNELING IN A 1D POTENTIAL' 
' ' 
' This file gives a graphic illustration of a particle tunneling ' 
' in a potential consisting of a barrier (or well) and a step.' 
' ' 
' It does so by utilizing the known analytic form of scattering  ' 
' solutions.' 
' ' 
' Summation over a finite number of Fourier frequencies (wave ' 
' numbers) allows us to follow the evolution of the state. '};

tt1=text0([.15 .5 .75 .35],txt);drawnow;

N=200; %  CHOOSE another value - depending on your computer.
nn=[-N:N];% 2*N+1 terms in Fourier series.
M=241;  % is the number of space points in [-5,7]; CHANGE HERE.
x=linspace(-5,7,M); 
%k0=300; % initial velocity
x0=4; % Initial position (average).
D=50; % D determines width of  wave packet, CHANGE HERE
delta=3*sqrt(D);
kk= delta*nn/N; % wave number lattice for incoming particle.
% Initial Gaussian wave packet
cn=exp(-kk.^2/D); 
nn0=sum(cn.^2); % Normalization.
% Fourier coeffs for initial wave packet-normalized = 
cn=cn.*exp(-i*x0*kk)/sqrt(nn0);
h=12/(M-1);

x1=-5:h:0-h; x11=ones(size(x1));
x2=0:h:2-h; x21=ones(size(x2));
x3=2:h:7; x31=ones(size(x3));

yy= exp(-i*kk'*x); % Each harmonic a row vector.
% Make a matrix where each each row is a harmonic function
% multiplied by the Fourier coefficient of the initial wave packet
yy1=zeros(size(yy'));
yy1(1,:)=cn;
yy1=cumsum(yy1);
yy1=yy1'.*yy; % Rows = space coord , columns = Fourier index
yy=abs(sum(yy1));
clear yy1
maxxy=max(yy);

% General form of the potential
y=[2*x11,10*x21,0*x31];
%Plot the general form of the potential and the initial wave packet:
cbutt(.75,.01);delete(tt1),
% pause
figure(gcf);


l1=plot(x,y,'r');
xlabel('Press the button to continue!');

hold on;
plot(x,yy);
title('This is the general potential and the initial wave packet (absolute value)');
set(l1,'LineWidth',2);
axis('off');
l1=text(3,5,'Incident wave packet');
set(l1,'FontName','palatino');
set(l1,'FontSize',18);
l1=text(-4.5,3,'V(x) = A');
set(l1,'FontName','palatino');
set(l1,'FontSize',18);
l1=text(0,11,'V(x) = B');
set(l1,'FontName','palatino');
set(l1,'FontSize',18);
l1=text(2,-.5,'x = 1');
set(l1,'FontName','palatino');
set(l1,'FontSize',18);
l1=text(-0.5,-.5,'x = 0');
set(l1,'FontName','palatino');
set(l1,'FontSize',18);

hold off

cbutt(.75,.01);
clear yy1 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q1=2; 

while  q1 >=  1

% INPUT 
if q1==2

str='1000';
tt2=text0([.15 .8 .75 .05],'Set the value of A');
gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume'),
AA=edit1([.75 .8 .15  .05],str);
A1=eval(AA); 
delete(tt2),


str='-1000';
tt2=text0([.15 .8 .75 .05],'Set the value of B');
gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume'),
AA=edit1([.75 .8 .15  .05],str);
A2=eval(AA);delete(tt2),


end 

str='50';
tt2=text0([.15 .8 .75 .05],'Set the velocity > 0');
gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume'),
AA=edit1([.75 .8 .15  .05],str);
k0=eval(AA);delete(tt2),
rbutt([.75 .01 .15 .05],'WAIT', ''),drawnow;


kk= k0 + delta*nn/N; % wave number lattice for incoming particle.

% Scale to get wave packet and the chosen potential into one picture:
max1=max([A1,A2,0]);
min1=min([A1,A2,0]);
range=max1-min1 + maxxy;
scale=0.6*maxxy/range;
offset=-min1/scale;
y=scale*([A1*x11,A2*x21,0*x31] -min1);

%% Calculate the scattering states which correspond to solutions in 
% [-infty, 0] which travel to the left or decay exponentially in that direction:

k1=ones(size(kk));
ee=kk.^2/2 +eps; % Kinetic energy
ee2=ee - A2; % ditto in [0,2]
ee1=ee - A1; % ditto in [-5,0]. 
ss2=sign(ee2);ss2p=(ss2+1);ss2n=(ss2-1);
ss1=sign(ee1);ss1p=(ss1+1);ss1n=(ss1-1);
kk2=sqrt(ss2p.*ee2)+eps - i*sqrt(ss2n.*ee2); % wave number in [0,1]
kk1=sqrt(ss1p.*ee1)+eps - i*sqrt(ss1n.*ee1); % wave number in [-infty,0]

aa2=0.5*(1+kk1./kk2); bb2=0.5*(1-kk1./kk2); % BC in x = 0
aa=0.5*exp(-i*2*kk).*((1+kk2./kk).*aa2.*exp(+i*2*kk2) ...
+ (1-kk2./kk).*bb2.*exp(-i*2*kk2));  % BC in x = 1
bb=0.5*exp(+i*2*kk).*((1-kk2./kk).*aa2.*exp(i*2*kk2) ...
+ (1+kk2./kk).*bb2.*exp(-i*2*kk2)); 

%% Now the whole solution, avoid making too many large matrices
%% Glue the pieces together.

yy = exp(-i*kk1'*x1)./(aa'*x11);
w = exp(-i*kk2'*x2).*((aa2./aa)'*x21) + exp(i*kk2'*x2).*((bb2./aa)'*x21);
yy =[yy,w];
w = exp(-i*kk'*x3) + exp(i*kk'*x3).*((bb./aa)'*x31);
yy=[yy,w];
xx=[x1,x2,x3];

clear w  aa bb aa2 bb2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make a matrix where each each row is a harmonic function 
% multiplied by the Fourier coefficient of the initial wave packet. 

yy=(cn'*ones(1,M)).*yy; % Rows = space coord, columns = Fourier index.


% Get the scales.
plot(xx,abs(sum(yy)));
axis('off'); xy=axis;


% We choose a sufficiently large unit of time to display the translation 
% of the wave packet.
tau=1/6/abs(k0);

rr=maxxy;
time=0;
figure(gcf);

		while rr > 0.1*maxxy

time=time+1;
yyt=abs(exp(-i*time*tau*kk.^2/2)*yy);% 
l1=plot(xx,[y;yyt]);set(l1,'LineWidth',2);
axis(xy),axis('off'),
text(0.8,0.8*xy(4),sprintf('time = %g',time));
drawnow
			if time < 100
			rr=max(yyt);
			else
			rr=0;
			end

		end 

clear yy % Get rid of the largest matrix!!!

rbutt([.15  .01  .15 .05],'New speed','uiresume, obutt, q1=1;'); 
rbutt([.3  .01   .15  .05],'New A,B','uiresume, obutt,  q1=2;');
bbutt([.45  .01  .15 .05],'BACK','close;q1=3;'); 
bbutt([.6  .01  .15  .05],'MAIN MENU','close;q1=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','uiresume;q1=0;'); 
uiwait;

		if q1==4
		start; return;

		elseif q1==3
		wavepac; return;

		end 

end  % end of q1,
close;

disp('> Type <tunnel> to do this all over again! ');



