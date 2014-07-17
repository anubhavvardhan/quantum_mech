%> <hard.m> calculates the partial cross sections for a hard sphere.
%> The program asks for a maximal angular momentum Lmax and displays 
%> the successive sums partial cross sections up to Lmax. Note that 
%> the convergence of the sum over L is very slow. Do not choose
%> too large a value of Lmax, say over 100.
%> The theory can be found in Messiah Ch X.13
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

clear;close; disp('> Welcome to <hard>');

q=1;

while  q==1

q1=1;

while q1==1

txt={' HARD SPHERE SCATTERING IN 3D' 
' ' 
' This file calculates the partial cross sections for a hard sphere.' 
' The program asks for a maximal angular momentum Lmax and displays ' 
' the successive sums of partial cross sections up to Lmax. ' 
' The convergence of the sum over L is slow. Do choose L < 100.' 
' ' 
' The program then displays the differential cross section as a' 
' function of the scattering angle theta as a function of momentum,' 
' both the asymptotic value and the finite sum over the parial' 
' cross sections up to Lmax.' 
' ' 
' Now set the value of Lmax '};

% figure(gcf),

tt1=text0([.15 .45 .75 .4 ], txt);
title('Hard sphere scattering'),axis('off'),
str='21';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .45 .15 .05], str);
Lmax=eval(ee); delete(tt1);

kmax=Lmax+5;
kr0= kmax*linspace(0.001,1,60); % wave no in multiples of x00 = k*r0
kr01=ones(size(kr0));
kr00= zeros(size(kr0));

w=[];

L=[0:Lmax]';% range of angular momenta
L1=ones(size(L));
 
% using standard Bessel function routine
norm=L1*sqrt(pi./2./kr0);
bess1=besselj(L+0.5,kr0)'.*norm;
bess2=bessely(L+0.5,kr0)'.*norm;

norm2=bess1.^2 + bess2.^2 +eps;

pwave=(bess2 - i*bess1).*bess1./norm2;
pwave=pwave./(L1*kr0);
phase=atan(bess1./bess2);
phase=glue(phase,pi);
plot(kr0,phase), axis tight,
title('Partial wave phase shifts for hard sphere');
xlabel('Wave number');

sigma0=4*pi*(2*L+1)*kr01.*(bess1./(L1* kr0)).^2 ./norm2;
sigma=[cumsum(sigma0);2*pi*kr01];

cbutt(.75,.01);

plot(kr0,sigma0), axis tight,
title('Partial cross sections for hard sphere'),
xlabel('Wave number'),cbutt(.75,.01)


plot(kr0,sigma), axis([0  kmax  0 4*pi]),
title('Summing the partial cross sections for hard sphere'),
xlabel('Wave number'), ylabel('Cross section in units of r0^2.'),
text(0.35*kmax, 3*pi,'The asymptotic cross section is 2 * pi * r0^2.');
text(0.3*kmax, 3.3*pi,'The curves are partial sums of the partial cross sections.');

rbutt([0.6  0.01  0.15  0.05],'New Lmax','uiresume,  q1=1;'),

gbutt([0.75  0.01   0.15  0.05],'CONTINUE','uiresume, q1=2;'),

uiwait;

obutt,

end
theta=linspace(eps,pi/Lmax,50)*pi;
zz=cos(theta);
theta1=ones(size(theta));
y =legpol(Lmax,zz);
LL=L*theta1;
yl=(2*LL+1).*y;
amp = pwave'*yl;
aamp=conj(amp).*amp;
maxa=max(max(aamp));

amp0= .25*(1+(cot(theta/2) .*besselj(1,kmax*sin(theta))).^2);

%% Plotting the cross sections
l1=plot(theta,amp0); axis('tight'),
set(l1,'LineWidth',2),
title('Asymptotic form of cross section'),
xlabel('Angle \theta '), ylabel('Cross section'),
legend('Asymptotic form');
hold on 
cbutt( .75, .01),
plot(theta,aamp),
title('Angular dependence of \sigma(\Omega) at different values of kr_{0}'),
hold off

cbutt( .75, .01),

%% Plotting the differential cross sections
l1=plot3(theta,kmax*kr01,amp0); hold on
set(l1,'LineWidth',2),
mesh(theta, kr0, aamp), axis('tight'),
xlabel('Angle \theta'), ylabel('Momentum kr0'),
view([1 -1 3]); hold off,

rbutt([0.3  0.01  0.15  0.05],'REPEAT','close,q=1;'); 
bbutt([0.45  0.01  0.15  0.05],'BACK','close,q=2;'); 
bbutt([0.6  0.01  0.15  0.05],'MAIN MENU','close,q=3;'); 
bbutt([0.75   0.01   0.15  0.05],'QUIT','close,q=4;'); 
uiwait; 

end

		if q==2

		scatt3d; return;
		
		elseif q==3
		
		start;return;		
		
		end

disp('> Type <hard> to start all over again!');




