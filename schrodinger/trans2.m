
%> The file <trans2.m> calculates the transmission coefficient 
%> for a finite series of quare potential barriers with a 
%> variable linear bias. 
%> It uses the method of multiplying transfer matrices.
%> The energy range can be chosen interactively.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear;close; disp('> Welcome to <trans2>!');

txt={' TRANSMISSION THROUGH A BIASED SET OF POTENTIAL BARRIERS' 
' ' 
' The transmission coefficient is calculated for a finite set of ' 
' square potential barriers, with a linear bias, as displayed in ' 
' the figure below. The width of each barrier is set = 1.' 
' The particles come in from the right' 
' Choose the number N and the height H of the barriers, the distance' 
' D between them, and the bias B.' 
''};

t1 = text0([.15  .5 .75 .35], txt);
subplot(2,1,1); axis('off');
subplot(2,1,2);axis('off');
l1=line([0,1],[0,0]);set(l1,'LineWidth',3);
l1=line([1,1],[0,1]);set(l1,'LineWidth',3);
l1=line([1,2],[1,1.05]);set(l1,'LineWidth',3);
l1=line([2,2],[.05,1.05]);set(l1,'LineWidth',3);
l1=line([2,5],[0.05,0.2]);set(l1,'LineWidth',3);
l1=line([5,5],[0.2,1.2]);set(l1,'LineWidth',3);
l1=line([5,6],[1.2,1.25]);set(l1,'LineWidth',3);
l1=line([6,6],[.25,1.25]);set(l1,'LineWidth',3);
l1=line([6,9],[.25,.4]);set(l1,'LineWidth',3);
l1=line([9,9],[.4,1.4]);set(l1,'LineWidth',3);
l1=line([9,10],[1.4,1.45]);set(l1,'LineWidth',3);
l1=line([10,10],[.45,1.45]);set(l1,'LineWidth',3);
l1=line([10,11],[.45,.45]);set(l1,'LineWidth',3);

string='3';
t2 = text0([.15  .55 .4 .05], 'Set the number N of barriers');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.75 .55 .15  .05],string);delete(t2),
N=eval(ee);

string='10';
t2 = text0([.15  .55 .4 .05], 'Set the height H of the barriers');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.75 .55 .15  .05],string);delete(t2),
H=eval(ee);

string='3';
t2 = text0([.15  .55 .4 .05], 'Set the distance D between barriers');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.75 .55 .15  .05],string);delete(t2),
D=eval(ee);

string='3';
t2 = text0([.15  .55 .4 .05], 'Set the total bias B');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.75 .55 .15  .05],string);delete(t2),
B=eval(ee);

P='You have chosen N = N0, H = H0, D = D0, B = B0';
P=sprintf(strrep(P,'N0','%g'),N); 
P=sprintf(strrep(P,'H0','%g'),H); 
P=sprintf(strrep(P,'D0','%g'),D); 
P=sprintf(strrep(P,'B0','%g'),B); 
t2 = text0([.15  .55 .6 .05], P);

cbutt(.75, .01);delete(t1,t2),

rbutt([.75 .01 .15 .05],'WAIT..',''),drawnow;
q0=2;

while q0==1|q0==2

if q0==2

X=[0,1]; V=H; 
%  V=0;
E0=[.1*V,5*(V+B)]; wav0=sqrt(2*E0 +eps);

q0=1;

end 

N00=300; %%% Change here
wave1=linspace(wav0(1), wav0(2),N00); %% incoming wave no

E=wave1.^2/2; E1=ones(size(E));

E=E+B;
G=B/(N + (N-1)*D);

wave2=sqrt(2*E); %% transmitted wave no
ww1=trf2([0,1],[H,H+G],E);
% ww0=trf1(D*X,0,E);
% ww2=mult4(ww0,ww1);ww=ww2;
ww=ww1;

 		for  iter=1:N-1

		L=1+(iter-1)*(D+1);
		ww1=trf2([0,1],H+G*[L,L+1],E);
		ww2=trf2([0,D],G*[L+1,L+D+1],E);
		ww1=mult4(ww1,ww2);
		ww=mult4(ww1,ww);

		end

f=ww(1,:)-i*wave2.*ww(3,:);
Df=ww(2,:)-i*wave2.*ww(4,:);

DD = f + i*Df./wave1;
% NN = f - i*Df./waveno;
TT=2./DD;% this is the transmitted amplitude
% RR=NN./DD;% this is the reflected amplitude
TC=abs(TT).^2.*wave2./wave1; % the transmission coefficient

wave1=wave1/sqrt(2*V);
subplot(2,1,1), 
xy=[wave1(1) wave1(N00)  0  1];
plot(wave1,TC); axis(xy);
title('Transmission coefficient as a function of normalized wave number');

z=.5*(ww(1,:) + ww(4,:)); z=crop(z,5,-5);

subplot(2,1,2)
xy=[wave1(1) wave1(N00) -5  5];
plot(wave1, z);axis(xy);hold on
plot(wave1,[E1;-E1],'r');hold off
title('Trace of transfer matrix (x 0.5), as a function of normalized wave number');

rbutt([.15  .01  .15  .05],'ZOOM IN','uiresume');
rbutt([.3   .01  .15  .05],'PANORAMA','uiresume, q0=2;');
rbutt([.45  .01  .15  .05],'FROM TOP','uiresume, q0=0;'); 
bbutt([.6   .01  .15  .05],'BACK','uiresume, q0=-1;'); 
bbutt([.75  .01  .15  .05],'QUIT','close, q0=-2;'); 
uiwait;

if q0==1
[x,irrelevant] = ginput(2); % gets 2 points from the current axes
wav0=x*sqrt(2*V);

elseif q0==0

trans2; return

elseif q0==-1

scatt; return

end

end

disp('> Type <trans2> to repeat this!');
