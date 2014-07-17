
%> The file <pw2d.m> treats 2D partial wave scattering, calculating
%> the scattering on a circular 'hard cylinder', and on an attractive
%> circular potential of arbitrary constant depth.
%> 
%> Choose radius R and potential depth D, and the energy range.
%> The partial cross sections and their sum are displayed. 
%> Reference:
%> SK Adhikari, American Journal of Physics, 54(4), 362 (1986)
%>
%> © Goran Lindblad - gli@theophys.kth.se

% solving for scattering amplitudes in 2D???
% D^2 u = 2(V-E)u
% solution are 
% J_m/Y_m (kr) cos/sin( m phi)   for E > V  k^2 = - 2(V-E)
% where J_m is regular at the origin
% I_m(kappa r) cos/sin(m phi) for E < V   kappa^2 = 2(V-E) and r < a
% regular at the origin
% Here we need the cos(m phi) compontents only
% DJ_m = - J_(m+1) + m J_m/z for J_m, Y_m as functions of z 

clear, close; disp('> Welcome to <pw2d>!');


txt={' CIRCULAR SQUARE WELL SCATTERING IN 2D' 
' ' 
' The script <pw2d> calculates the 2D scatterin on a square well potential' 
' of cylindrical symmetry and arbitrary depth.' 
' ' 
' V(r) =  - D  < 0  for 0 < r < R, V = 0 for r > R. ' 
' ' 
' The partial cross sections and their sum will be displayed. ' 
' You have to choose R, D and the wave number range.' };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED PARAMETERS, CHANGE HERE
Nk=500; %% no of lattice points in wave number
m_max=40; %% maximal angular momentum
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tt1=text0([.15 .5 .75 .35], txt);

tt2=text0([.15 .5 .75 .05], 'Choose the pair [R, D] >>> ');

str='[3,1.5]';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .5 .15  .05],str);
ee=eval(ee); delete(tt2), 
R=ee(1); D=ee(2);

text2='Choose the wave number range [kmin,kmax] >>> ';
tt2=text0([.15 .5 .75 .05],text2 );

str='[0,1.5]';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .5 .15  .05],str);
ee=eval(ee); delete(tt1,tt2), 
kmin=max(ee(1), .01); %% nonzero min to avoid divergences
kmax=ee(2);


kk=linspace(kmin,kmax,Nk); %% lattice of wave numbers
k1=sqrt(kk.^2 +D); %% wave number in potential V = - D;

mm=[0:m_max]; %% angular momenta


theta= linspace(0,1)*pi;
rad=linspace(1,30,50);

bess=besselj(mm, rad');
coss=cos(mm'*theta);
coeff=2*exp(i*mm*pi/2);
coeff(1)=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now solve for the phase shifts

% Calculating logarithmic derivative at r = R using
% DJ_m = - J_(m+1) + m J_m/z for J_m, Y_m as functions of z 

% Solution inside, only regular one!

zz=k1*R;
[mmm, kk1]=meshgrid(mm,k1);

bj=besselj([0,mm+1],zz');
bj1 = bj(:, mm+1); 
bj2 = bj(:, mm+2);

% LD= mmm/R - bj1./bj2; 
%% Then solutions outside

zz=kk*R;

bj=besselj([0,mm+1],zz');
by=bessely([0,mm+1],zz');

[mmm, kkk]=meshgrid(mm,kk);
[mmm,zzz]=meshgrid(mm,zz);

nn=kkk.^mmm; %% normalization factor to avoid divergences

cj1 = bj(:, mm+1)./nn;
cj2 = bj(:, mm+2)./nn;
cy1 = by(:, mm+1).*nn; 
cy2 = by(:, mm+2).*nn;

%% Now the CC at r = R gives the formula for phase shift
%% tand = tan(delta) as a function of m and k

nom=kk1.*bj2.*cj1- bj1.*cj2.*kkk;
dom=kk1.*bj2.*cy1 - bj1.*cy2.*kkk;
tand =  nn.^2.*nom./dom;

%% SIGN?????????

tand=tand'; %% transpose

del=atan(tand);
del=glue(del,pi);
plot(kk',del),
title('Phase shifts for m = 0, 1, 2 ,...');
cbutt(.75,.01), 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now calculate partial and total cross sections

sin2d = tand.^2./(1 +tand.^2);

ee=2*ones(size(mm));
ee(1) = 1;
sig= (ee'*(4./kk)).*sin2d;

sigtot=sum(sig);
sigmax=max(max((sigtot)));

plot(kk,sigtot,'k'), 
hold on,
plot(kk,sig), axis tight
title('Partial and total cross section versus wave number');
cbutt(.75,.01),
hold off;

n0 = round(Nk/5);
%% calculate scattering amplitude as function of theta
theta=linspace(0,1,200)*2*pi;

q1=1;

while q1==1

k0= kmin + n0*(kmax-kmin)/Nk;
del1=del(:,n0);

amp= exp(i* del1).*sin(del1);
costheta = cos(theta'*mm);

amp=sqrt(2/pi)*costheta*diag(ee)*amp;
% plot(theta',abs(amp));
f1=polar(theta'+pi/2,abs(amp));
set(f1,'LineWidth',3)
%axis tight;
title('Polar diagram of abolute value of scattering amplitude');
xlabel('Forward direction at top. Change wave number with slider!');

gbutt([0.75  0.01  0.15  0.05],'PASS','q1=0;uiresume');
g0 =slider0([.92  .12  .05  .8],[0,1,n0/Nk],'q1=1; uiresume;');
uiwait; x=get(g0,'Value');   delete(g0);
n0=max(ceil(x*Nk),1);

end


theta=linspace(0,1,200)*pi;
%% Scattered wave, 
rr=linspace(R, 5*R, 100);
[mmm,rrr]=meshgrid(mm',rr);

jv=besselj(mm',k0*rr);
yv=bessely(mm',k0*rr);

xx=cos(theta')*rr;
yy=sin(theta')*rr;

% first display the incoming plane wave

% amp=cos(theta'*mm)*diag(ee.*exp(i*mm*pi/2))*jv';
% ramp=real(amp);
% surf(xx,yy,2*ramp); axis equal,
% view(2),
% cbutt(.75,.01)

q1=1;

subplot(2,1,1)
plot(kk,sigtot,'k'), 
hold on,
plot(kk,sig), axis tight
title('Partial and total cross section versus wave number');


while  q1==1

subplot(2,1,1);

p2=plot(kmin+n0*(kmax-kmin)/Nk,sigmax*.9,'^ r');
set(p2,'LineWidth',3);
xlabel('Choose wave number with slider!');
k0= kmin + n0*(kmax-kmin)/Nk;
del1=del(:,n0);
jv=besselj(mm',k0*rr);
yv=bessely(mm',k0*rr);

% then display the scattered amplitude

coeffc=diag(ee.*exp(i*mm*pi/2).*sin(del1').*exp(i*del1'));
coeffs=diag(ee.*exp(i*mm*pi/2).*sin(del1').^2);
amp=cos(theta'*mm)*coeffc*(jv' + i*yv'); %%% sign?

ramp=abs(amp);
mamp=max(max(ramp));
ramp=R*ramp/mamp;


subplot(2,1,2)
mesh(xx,yy,ramp),axis equal,
%view([0, 1, 5]);
title('The absolute value of scattered wave for chosen wave number');
xlabel('The x-axis is the forward direction!');
view(2), rotate3d,% colormap('copper'),
gbutt([0.75  0.01  0.15  0.05],'PASS','q1=0;uiresume'); 
% cbutt(.75,.01),


g0 =slider0([.92  .12  .05  .8],[0,1,n0/Nk],'q1=1; uiresume;');
uiwait; x=get(g0,'Value');   delete(g0,p2);
n0=max(ceil(x*Nk),1);

end

txt='> Type <pw2d> to start all over again!';

rbutt([.3  .01  .15  .05],'REPEAT','uiresume, pw2d;return;'); 
bbutt([.45  .01  .15  .05],'BACK','close,scatt3d; return,'); 
bbutt([.6   .01  .15  .05],'MAIN MENU','close,start;return,'); 
bbutt([.75  .01  .15  .05],'QUIT','close;disp(txt)'); 
uiwait;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
