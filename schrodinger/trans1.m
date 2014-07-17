
%> The file <trans1.m> calculates the transmission coefficient 
%> for a finite series of quare potential barriers.
%> It uses the method of multiplying transfer matrices.
%> The energy range can be chosen interactively.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear, close; disp('> Welcome to <trans1>!');

txt={' TRANSMISSION THROUGH A FINITE SEQUENCE OF POTENTIAL BARRIERS' 
' ' 
' The transmission coefficient is calculated for a finite set of ' 
' square potential barriers, as displayed in the figure below.' 
' The width of each barrier is set = 1.' 
' ' 
' Choose the number and the height of the barriers, and the distance' 
' between them.' 
''};


subplot(2,1,1), axis('off'),
t1 = text0([.12  .5 .78 .35], txt);
subplot(2,1,2), axis('off'),
l1=line([0,0],[0,1]);set(l1,'LineWidth',3),
l1=line([0,1],[1,1]);set(l1,'LineWidth',3),
l1=line([1,1],[0,1]);set(l1,'LineWidth',3),
l1=line([1,4],[0,0]);set(l1,'LineWidth',3),
l1=line([4,4],[0,1]);set(l1,'LineWidth',3),
l1=line([4,5],[1,1]);set(l1,'LineWidth',3),
l1=line([5,5],[0,1]);set(l1,'LineWidth',3),
l1=line([5,8],[0,0]);set(l1,'LineWidth',3),
l1=line([8,8],[0,1]);set(l1,'LineWidth',3),
l1=line([8,9],[1,1]);set(l1,'LineWidth',3),
l1=line([9,9],[0,1]);set(l1,'LineWidth',3),
l1=line([9,11],[0,0]);set(l1,'LineWidth',3);

string='3';
t2 = text0([.15  .55 .4 .05], 'Set the number N of barriers');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.7 .55 .2  .05],string);delete(t2),
N=eval(ee);

string='10';
t2 = text0([.15  .55 .4 .05], 'Set the height H of the barriers');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.7 .55 .2  .05],string);delete(t2),
H=eval(ee);

string='3';
t2 = text0([.15  .55 .4 .05], 'Set the distance D between barriers >>> ');
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
ee=edit1([.7 .55 .2  .05],string);delete(t2),
D=eval(ee);

P='You have chosen N = N0, H = H0, D = D0';
P=sprintf(strrep(P,'N0','%g'),N),
P=sprintf(strrep(P,'H0','%g'),H), 
P=sprintf(strrep(P,'D0','%g'),D), 
t2 = text0([.15  .55 .5 .05], P);

cbutt(.75, .01);delete(t1,t2),

q0=2;

while q0==1|q0==2

if q0==2

X=[0,1]; V=H; 
%  V=0;
E0=[.3*V,5*V]; wav0=sqrt(2*E0 +eps);

q0=1;

end 

N0=300; %%% Change here!
waveno=linspace(wav0(1), wav0(2),N0);
E=waveno.^2/2;E1=ones(size(E));

ww1=trf1(X,V,E);ww0=trf1(D*X,0,E);
ww2=mult4(ww0,ww1);ww=ww2;

 		for  iter=1:N-1

		ww=mult4(ww2,ww);
		
		end

f=ww(1,:)-i*waveno.*ww(3,:);
Df=ww(2,:)-i*waveno.*ww(4,:);

DD = f + i*Df./waveno;
% NN = f - i*Df./waveno;

TT=2./DD;% this is the transmitted amplitude
% RR=NN./DD;% this is the reflected amplitude
TC=abs(TT).^2; % the transmission coefficient
% close(f1)

waveno=waveno/sqrt(2*V);
xy=[waveno(1) waveno(N0) 0  1];

subplot(2,1,1)
plot(waveno,TC);axis(xy);
title('Transmission coefficient as a function of normalized wave number');
z=.5*(ww(1,:) + ww(4,:));
xy=[waveno(1) waveno(N0) -5  5];
subplot(2,1,2)
plot(waveno, z);axis(xy); hold on;
plot(waveno,[E1;-E1],'r');hold off
title('Trace of transfer matrix (x 0.5), as a function of normalized wave number');

rbutt([.15  .01  .15  .05],'ZOOM IN','uiresume');
rbutt([.3  .01  .15  .05],'PANORAMA','uiresume, q0=2;');
rbutt([.45  .01 .15  .05],'FROM TOP','uiresume, q0=0;'); 
bbutt([.60  .01  .15  .05],'BACK','uiresume, q0=-1;'); 
bbutt([.75  .01  .15  .05],'QUIT','close, q0=-2;'); 

uiwait;

if q0==1
[x,irrelevant] = ginput(2); % gets 2 points from the current axes
wav0=x*sqrt(2*V);

elseif q0==0

trans1; return

elseif q0==-1

scatt; return

end

end

disp('> Type <trans1> to repeat this!');
