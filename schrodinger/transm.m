
%> The file <transm.m> calculates the transmission coefficient in 1D 
%> potential scattering using the Numerov algorithm for integration. 
%> The scattering potential is nonzero only on [0,1]. 
%> A lattice of values for the energy E is chosen in the CW.
%> The particle comes in from x = + infinity, we start the integration
%> from the transmitted wave in x = 0, wave number sqrt(2*E). 
%> At x = 1 we adapt the incoming and reflected waves and renormalize 
%> to find the transmission coefficient, which is displayed as a 
%> function of wave number. Forbidden wave number bands are displayed
%> for a choice of unit cell spacings.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125

clear; close; disp('> Welcome to <transm>');
 
txt={ ' TRANSMISSION COEFFICIENTS IN 1D SCATTERING	' 
' ' 
' This file calculates the energy-dependence of the transmission ' 
' coefficient for the 1D Schrodinger equation using the Numerov ' 
' algorithm.' 
' ' 
' The scattering potential is chosen on the next page.' 
' ' 
' The transmission coefficient is displayed as a function of' 
' wave number. After a choice of lattice spacing the allowed' 
' and forbidden energy bands are displayed.'};

tt1=text0([.15  .45 .75 .4], txt);

cbutt(.75 , .01), delete(tt1);

Q1=4; 

	while Q1 > 0
   
			if Q1==4

choice;  %% calls <choice.m> for choice of potential

XX=linspace(0,1); 
subplot(2,2,3)
p1=plot(XX,pot1(XX)), axis('square'),
set(p1,'LineWidth',4),
title('This is your potential, normalized.'),
xlabel('Press the button to continue!');

txt={'You must choose the number of points in the space grid,'
'the number of points in the energy grid, the multiplier of the'
'the potential (including sign!), and the energy interval..'};
tt1=text0([.15 .5 .75 .15], txt);

str='500';
tt2=text0([.15 .5 .75 .05],'Give the number of space grid points ');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .5 .15  .05],str);
N=eval(ee);delete(tt2);

str='200';
tt2=text0([.15 .5 .75 .05],'Give the number of energy grid points ');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .5 .15  .05],str);
M=eval(ee);delete(tt1,tt2);

			Q1=3; 		
			
			end

			if Q1==3
			
str='-300';
tt2=text0([.15 .5 .75 .05],'Give the multiplier of the potential (with sign) ');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .5 .15  .05],str);
multiplier=eval(ee);delete(tt2);
			
subplot(2,2,3),
p2= plot(XX,multiplier*pot1(XX));
set(p2,'LineWidth',4);
title('This is your potential, with multiplier')

%%%%%%%%%%%%%%%%%%%%%%

h=1/(N-1); X=[0:h:1]'; uniX=ones(size(X));
V=multiplier*pot1(X);

			Q1=2;		

			end

			if Q1==2			
			
%%%%%%%%%%% INPUT %%%%%%%%%%%

str='[0,300]';
tt2=text0([.15 .5 .75 .05],'Choose the energy interval');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .5 .15  .05],str);
energy=eval(ee);delete(tt2);
waveno=linspace(sqrt(2*(energy(1)+eps)),sqrt(2*energy(2)),M);
Y=0.5*waveno.^2; % Energy
Y1= ones(size(Y)); Y0= zeros(size(Y));
X1=Y1;	X2=-waveno*i;

%%%%%%% NUMEROV INTEGRATION %%%%%%%%%%%%%
%ww=numerov1(0,1,N,'pot1(X)',multiplier,Y,X1,X2);

ww=transfer(0,1,N,'pot1(X)',multiplier,Y);

f=ww(1,:)-i*waveno.*ww(3,:);
Df=ww(2,:)-i*waveno.*ww(4,:);
DD = f + i*Df./waveno;
NN = f - i*Df./waveno;

TT=2./DD;% this is the transmitted amplitude
RR=NN./DD;% this is the reflected amplitude
TC=abs(TT).^2; % the transmission coefficient

subplot(1,1,1); plot(waveno,[TC ; Y1-TC - abs(RR).^2]),
axis('tight'), 
title('Transmission coefficient as a function of wave number');
xlabel('Wave number. - Press button to continue!');

rbutt([0.55 ,0.01,0.2,0.05],'Bloch solutions','uiresume, Q1=1;'); 
gbutt([0.75 ,0.01,0.15,0.05],'CONTINUE','uiresume; Q1=2;');
uiwait;
			
			end 

			if Q1==1
	
str='2';
tt2=text0([.15 .5 .75 .05],'Choose the unit cell legth');
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .5 .15  .05],str);
L=eval(ee);delete(tt2);		
%% Trace of the transfer matrix for a unit cell:
MM=0.5*((ww(1,:)+ww(4,:)).*cos(waveno*(L-1)) +...
 ((ww(2,:)./waveno) - ww(3,:).*waveno).*sin(waveno*(L-1)));


MM=crop(MM,1,-1); %% Restricts the values of MM to [-1,1].
y=acos(MM); %% Bloch wave number.
w2=[-0.1,pi,min(waveno),max(waveno)];
MM0 = 1-abs(MM); %% This is zero where MM has values +/- 1.
MM0=~MM0; %% = 1 where MM has values +/- 1, else zero.
% % MM0(M)=1;MM0(1)=1;

plot(y,waveno);axis(w2);
xlabel('Bloch wave number K*L');
ylabel('Wave number of Bloch solution and allowed bands');
title('Wave number of Bloch solutions versus Bloch wave number');
hold on;  barh(waveno,-0.1*MM0); hold off;

end

%%%%%%%%%%% INPUT %%%%%%%%%%%
rbutt([0.35  0.01   0.2  0.05],'New unit cell','uiresume; Q1=1;'); 
rbutt([0.55  0.01  0.2  0.05],'- energy range','uiresume; Q1=2;');
gbutt([0.75  0.01   0.15  0.05],'CONTINUE','uiresume; Q1=0;');
uiwait;

if Q1~=0
obutt;
else

rbutt([0.2  0.01   0.2  0.05],'New multiplier','uiresume; Q1=3;');
rbutt([0.4  0.01  0.2  0.05],'- potential','close; Q1=4;');
bbutt([0.6  0.01   0.15  0.05],'BACK','close;Q1=5;'); 
bbutt([0.75  0.01  0.15 0.05],'QUIT','close;Q1=0;'); 
uiwait;

if Q1==3
obutt;
end

end 
	if Q1==5
	
	per1d; return;

	end;
	
	end %%of Q1	
	
disp('> Type <transm> to do this problem again. Bye!');


