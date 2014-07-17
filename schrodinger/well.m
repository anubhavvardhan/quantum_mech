
%> The file <well.m> calculates the partial cross sections for a spherical
%> square well potential of arbitrary strength, including arbitrary sign. 
%> V(r) =  V for  r < A, V = 0 for r > A.
%> Choose A and V and the energy range.
%> The partial cross sections and their sum are displayed as 
%> functions of wave number, then the differential cross section 
%> as a function of the scattering angle for a chosen energy 
%> Reference: Messiah Ch IX--X.
%>
%> © Goran Lindblad - gli@theophys.kth.se


% GL 961125 

clear, close; disp('> Welcome to <well>!');

q=1;

while q==1

txt={' SPHERICAL SQUARE WELL SCATTERING' 
' ' 
' This file calculates the partial cross section for a square well' 
' potential of spherical symmetry and arbitrary strength and sign.' 
' ' 
' V(r) =  V for 0 < r < A, V = 0 for r > A. ' 
' ' 
' The partial cross sections and their sum are displayed. ' 
' You have to choose V, A and the wave number range.' };

tt1=text0([.15 .5 .75 .35], txt);

tt2=text0([.15 .5 .75 .05], 'Choose V, with sign!');

subplot(2,2,3);
axis([-.5 3.1 -3 2]);axis('off');axis('equal')
l1=line([0,3],[0,0]);
set(l1,'LineWidth',1);
l1=line([0,0],[1.3,-2]);
set(l1,'LineWidth',3);
l1=line([0,2],[-2,-2]);
set(l1,'LineWidth',3);
l1=line([2,2],[-2,0]);
set(l1,'LineWidth',3);
l1=line([2,3],[0,0]);
set(l1,'LineWidth',3);
t1=text(2.7,.2 ,'r');
set(t1,'FontSize',18);
t1=text(2,.2 ,'A');
set(t1,'FontSize',18);
t1=text(1,-1.7 ,'V');
set(t1,'FontSize',18);

str='-2';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .5 .15  .05],str);
V=eval(ee); delete(tt2), 

tt2=text0([.15 .5 .75 .05], 'Choose the radius A > 0');
str='2';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .5 .15  .05],str);
A=eval(ee); delete(tt1,tt2),

		if A < 0 

		disp('> The radius must be positive, changing sign!');	A=-A;

		else

		end

q=2; 

while  q == 2

tt2=text0([.15 .5 .75 .05], 'Choose the wave number range [k0,k1]');
str='[0,2]';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .5 .15  .05],str);
kkk=eval(ee); delete(tt2),

	if V > 0

Lmax=ceil(A*kkk(2)+2);    %  ceil(sqrt(2*E(2))); %choice of maximal angular momentum
%kmax = Lmax +2 %choice of maximal wave number

	else
Lmax=ceil(A*sqrt(kkk(2)^2 - 2*V)+2);  %choice of maximal angular momentum
%kmax=3*floor(sqrt(Lmax)) + 5; %choice of maximal wave number
	end

% disp(sprintf('> The maximal angular momentum is Lmax =%g ',Lmax));

L=[0:Lmax]'; %the range of angular momenta
oneL=ones(size(L));

% kkk=sqrt(2*E);
%%%%first the free space region
N=200; % The number of points in wave number range, change here!
epsilon=.0001; % cutoff for small wave numbers, change here!
kk0=linspace(kkk(1),kkk(2),N)+epsilon;
% the wave number range for the free particle
onek=ones(size(kk0));
bb01=besselj(L+0.5,kk0*A)';%.*norm0;
bb02=bessely(L+0.5,kk0*A)';%.*norm0;
bb03=besselj(L+1.5,kk0*A)';%.*norm0;
bb04=bessely(L+1.5,kk0*A)';%.*norm0;

% logarithmic derivatives
% logd01 =L*onek - (oneL*kk0).*bb03./bb01;
% logd02 =L*onek - (oneL*kk0).*bb04./bb02;
% logd01 = (L+1)*onek - A*(oneL*kk0).*bb03./bb01;
% logd02 = (L+1)*onek - A*(oneL*kk0).*bb04./bb02;
%%%%%%%% CHECK THIS FORMULA

logd01 = ((L+1)*onek).*bb01 - A*(oneL*kk0).*bb03; %non-normalized !!
logd02 = ((L+1)*onek).*bb02 - A*(oneL*kk0).*bb04;

%%%%%%now the potential well

kk=kk0.^2 - 2*V; % changes sign if V > 0
sigg=sign(kk);
sig2=sign(sigg+1);% 1's where kk is non-negative, zeros elsewhere.
sig1=onek-sig2; % 1's where kk is negative, zeros elsewhere.

kk1=sig1.*sqrt(-sig1.*kk+eps); % low wave nos  where kk is negative
kk2=sig2.*sqrt(sig2.*kk +  eps); % high - " -  positive

bb11=(oneL*sig1).*besseli(L+0.5,kk1*A)'+ (oneL*sig2).* besselj(L+0.5,kk2*A)';%.*norm1;
bb12=(oneL*sig1).*besseli(L+1.5,kk1*A)' + (oneL*sig2).* besselj(L+1.5,kk2*A)';%.*norm1;

%calculate the logarithmic derivatives

logd1 = ((L+1)*onek - A*oneL*(kk1+kk2).*(oneL*sigg).*bb12./bb11);

%the boundary condition in r=1 %% r=A 
tandelta=(logd01 - bb01.*logd1)./(logd02 - bb02.*logd1);
delta=atan(tandelta);

% We plot the phase shifts 
delta=glue(delta,pi);
subplot(2,1,1),f0=plot(kk0,delta);axis tight; hold on;
f1=plot(kk0,.5*pi*[1;-1]*ones(size(kk0)), 'k'); hold off;
title('The partial wave phase shifts as a function of wave number');

%calculate sin^2 delta and the cross section
sindelta2= tandelta.^2./(1+ tandelta.^2);
sigma = 4*pi*((2*L+1)*(onek./kk0./kk0)).*sindelta2;

sigsum=sum(sigma');
sigsumsum=sum(sigsum);
weight=sigsum/sigsumsum;
weight=[0:length(weight)-1; weight]';
% disp('> The relative contributions of the partial waves are:');
% disp(weight);

string=sprintf('Partial cross sections and their sum up to Lmax = %g ',Lmax);

% plot the lot, first the sum of the partial cross sections
cross=sum(sigma); mm1=max(cross);
subplot(2,1,2),
l0=plot(kk0,[cross;zeros(size(kk0))],'k');
axis tight;
set(l0,'Linewidth',2);
title(string);
xlabel('Wave number.');

%now plot the partial cross sections on the same scale
hold on 
plot(kk0,sigma);

cbutt(.75,.01);

string2=str2mat('> Now we will calculate the differential cross section',...
'> using the same parameter V.');

ampl= ((2*L+1)*(1./kk0)).*(tandelta + i*tandelta.^2)./(1 + tandelta.^2);
theta=linspace(0,1,50)*pi;
pol=legpol(Lmax,cos(theta));
ampl=pol'*ampl;
difcros=real(conj(ampl).*ampl);
[row,col]=size(difcros);
mm2=max(max(difcros));

subplot(2,1,1);
p1=plot(theta, difcros(:,1)); axis([0 pi 0 mm2]),
title('Differential cross section, chose energy with slider !');
str1='Wave number =';
str2='x';
str=sprintf(strrep(str2,'x','%g'),kkk(1));
tt1=text0([.5 .5 .2 .04], str1); 
tt2=text0([.7 .5 .1 .04], str);
par=1;
hold on
subplot(2,1,2);
p2=plot(kkk(1)+par*kkk(2)/col,0,'^ r');
set(p2,'LineWidth',3);
% nn=round(col/30);

q1=1; par=1;

while q1 == 1
delete(p1,p2,tt2);
subplot(2,1,1)

str=sprintf(strrep(str2,'x','%g'),kkk(1)+par*kkk(2)/col);
p1=plot(theta, difcros(:,par));set(p1,'LineWidth',3);
tt2=text0([.7 .5 .1 .04], str);
subplot(2,1,2);
p2=plot(kkk(1)+par*kkk(2)/col,mm1*.9,'^ r');
set(p2,'LineWidth',3);
gbutt([0.75  0.01  0.15  0.05],'PASS','q1=0;uiresume'); 
g0 =slider0([.92  .12  .05  .8],[0,1,par/col],'q1=1; uiresume;');
uiwait; x=get(g0,'Value')+eps;   delete(g0);
par=round(x*col)+1;
par=min(col,par);

end

delete(p1,p2),delete(tt2),
% str=sprintf(strrep(str2,'x','%g'),kkk(2));
% subplot(2,1,1)
% p1=plot(theta, difcros(:,col));
% tt2=text0([.7 .5 .1 .04], str);
% set(p1,'LineWidth',3);
 hold  off

% cbutt(.75,.01); 
delete(tt1);

%  make a log scaling of the cross section for the 3D graphics
difcros=log(difcros+1); 
subplot(1,1,1);
mesh(theta,fliplr(kk0),rot90(difcros));
axis tight;  view([1.3,1,5]);
title('Differential cross section as a function of angle and wave number, logarithmic scale!');
ylabel('Wave number'); xlabel('Angle theta in [0,pi]');
%zlabel('Logarithmic scale!');

rbutt([.15  .01  .15  .05],'REPEAT','uiresume, obutt, q=1;'); 
rbutt([.3   .01  .15  .05],'Set [k0,k1]','uiresume, obutt, q=2;');
bbutt([.45  .01  .15  .05],'BACK','uiresume, q=3;'); 
bbutt([.6   .01  .15  .05],'MAIN MENU','uiresume, q=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q=0;'); 
uiwait;

	end % q=2 returns

	end  % q=1 returns

	if q==3

	scatt3d;	return;

	elseif q==4

	start;return;	

	end

disp('> Type <well> to start all over again!');


