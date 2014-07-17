
%> The file <barrier.m> calculates the energy dependent transmission 
%> coefficient for a 'square' barrier (of either sign!) in 1D. 
%> The coefficient is first calculated using the Numerov algorithm.
%> This result is then compared with the known analytical formula given 
%> e.g. in Messiah Ch III. The transmission coefficient is displayed as 
%> a function of the energy renormalized by the barrier height.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125

clear; close;

disp('> Welcome to <barrier>');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% FIXED PARAMETERS, YOU CAN CHANGE HERE
%%%%%%
N = 300; %% Number of steps in numerical integration
M = 200; %% Number of points in energy lattice
A = 6; %% Defines maximal energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

txt={' SCATTERING FOR A SQUARE POTENTIAL BARRIER' 
' ' 
' The file <barrier.m> calculates the energy dependent transmission ' 
' coefficient for a square barrier (of either sign!). It is first ' 
' calculated by integrating the 1D Schrodinger equation using the ' 
' Numerov algorithm.' 
' ' 
' This result is then compared with the known analytical formula ' 
' given in many QM textbooks, e.g. Messiah Ch III. ' 
' The transmission coefficient is displayed as a function of the' 
' energy renormalized by the barrier height.' 
' ' 
'  You have to choose the barrier height!'};

t1=text0([.15 .5 .75 .35], txt);
cbutt(.75 ,.01), delete(t1);

%%%%%% THE POTENTIAL %%%%%%%    
x=linspace(-0.5,1.5,200);
y=0.5*(sign(x) + sign(1-x));
PP='0.5*(sign(X) + sign(1-X))';%potential defined as  a string.

Q1=1;

V0='100';

l1=plot(x,y);set(l1,'LineWidth',3);axis('equal'), 
title('Shape of the potential barrier'),

txt={' Insert the multiplier, then make a return!'
	' You can choose it to be of either sign!'};
t2=text0([.15 .8 .75 .08],txt);

gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.75 .8 .15 .05],V0); delete(t2); V=eval(ee);

while Q1==1

V1=abs(V)+eps;

if V > 0

E0=max(eps, 1-10/(V+ eps)); energy=[E0,A]*(V);

else 

E0 = 0; energy=[E0,A]*abs(V);

end

h=1/(N-1);X=[0:h:1]'; uniX=ones(size(X));

waveno=linspace(sqrt(2*(energy(1)+eps)),sqrt(2*energy(2)),M);

wav0=waveno./sqrt(2*V1); Y=0.5*waveno.^2; % Kinetic energy

Y1= ones(size(Y));Y0= zeros(size(Y)); X1=Y1; X2=-waveno*i;

%%%%%%% NUMEROV INTEGRATION %%%%%%%%%%%%%

delta=0.0; %%% 
ww=numerov1(0-delta,1+delta,N,PP,V,Y,X1,X2);

DD = ww(1,:) + i*ww(2,:)./waveno;
NN = ww(1,:) - i*ww(2,:)./waveno;

TT=2./DD;		% this is the transmitted amplitude
RR=NN./DD;	% this is the reflected amplitude
T0=abs(TT).^2; 

%%%%%%%% ANALYTICAL FORMULA %%%%%%%

Z=2*(Y-V); 
Z1=floor(0.5*(sign(Z)+1)); % Z > 0
Z2=ceil(0.5*(1-sign(Z))); % Z < 0
SS = sin(sqrt(Z1.*Z)) + sinh(sqrt(-Z2.*Z));
DD=1+V^2*SS.^2./Y./(abs(Z)+eps)/2; % See for instance Messiah Ch III.
T1=1./DD;

%%%%%% PLOTTING %%%%%%%%

xy=[sqrt(E0)-.02   sqrt(A)  -.1  1];

l1=plot(wav0,[T0 ;T1;T0-T1]); set(l1,'LineWidth',2), axis(xy),
title('Absolute value of transmission coefficient as a function of wave number'),
xlabel('Wave number, using a scale renormalized to the barrier height = 1!'),

P=' The multiplier = V';
P=sprintf(strrep(P,'V','%g'),V);  %%% Replace string by value 
text0([.4  .25  .3  .05], P);
t1=legend('Numerical','Exact', 'N-E',0);

cbutt(.75 ,.01),
% delete(t1);

plot(Y/V1,1-T1./T0), title('Displays 1-Texact/Tnumeric'),

rbutt([.45  .01  .15  .05],'REPEAT','uiresume; Q1=1;'),
bbutt([.6  .01   .15  .05],'BACK','close; Q1=2;'),
bbutt([.75  .01 .15  .05],'QUIT','close; Q1=3;'), 
uiwait;

	if Q1==1

obutt,
txt={' Insert the new multiplier at right >'};
t2=text0([.15 .85 .75 .05],txt);
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
V0=edit1([.75 .85 .15 .05],V0);
delete(t2); V=eval(V0);

	elseif Q1==2

	scatt; return;

	end 

end

disp('> Type <barrier> to do this problem again. Bye!');



%%% © Göran Lindblad 1996
