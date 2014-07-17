
%> The file <diffrac.m> calculates interference patterns from a multiple slit.
%> It uses <fresn.m> to calculate the Fresnel integrals.
%> Input: number of slits, width of slits, distance between centers
%> renormalized wave length sqrt(d*wavelength), d = distance to screen
%> Refs: GÄHLER, R & ZEILINGER, A [GAHLER]
%> Wave-optical experiments with very cold neutrons
%> Am J Phys 59(4), 316--324 (1991)
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear;close; disp('> Welcome to <diffrac>!');

q0=1;

N=500; % No of points on screen, change here!

% figure(gcf),

txt={' DIFFRACTION FROM MULTIPLE SLITS' 
' ' 
' The file <diffrac.m> calculates interference patterns for a multiple' 
' slit experiment. It uses the function file <fresn.m> to calculate the' 
' Fresnel patterns.' 
' ' 
' Input: number of slits, width of slits, distance between centers, ' 
' renormalized wavelength =  sqrt(d*wavelength), ' 
' where d = distance to screen.' 
' You can zoom the screen using the mouse, and change the wavelength' 
' using the slider! '};

tt1=text0([.15 .5 .75 .35],txt);

subplot(2,2,3),
c=8;
y0=linspace(-c,c,N);axis tight,
f0=fresn(y0);
plot(real(f0),imag(f0)), 
title('Cornu spiral = Fresnel integrals C(x)+iS(x)'),
xlabel('Value of Fresnel integral C(x)'),
ylabel('Value of Fresnel integral S(x)')
axis equal,


subplot(2,2,4),
y0=linspace(-c,c/2,N);
f0=fresn(y0)-0.5-0.5*i; f0=conj(f0).*f0/2;
plot(y0,f0),title('Fresnel diffraction by a straight edge'),
axis tight,
xlabel('Distance from straight edge'),
cbutt(.75,.01); delete(tt1),

while q0==1

tt2=text0([.15 .55 .75 .05], 'Set the number of slits:' );
str='6';
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume'),
ee=edit1([.75 .55 .15  .05],str),
n0=eval(ee); delete(tt2),

tt2=text0([.15 .55 .75 .05], 'Set the width of each slit:' );
str='1';
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .55 .15  .05],str);
a=eval(ee); delete(tt2);

tt2=text0([.15 .55 .75 .05], 'Set the distance between the centers of the slits:'  );
str='8';
gbutt([.75 .01 .15 .05],'CONTINUE','uiresume');
ee=edit1([.75 .55 .15  .05],str);
b=eval(ee); delete(tt2);


q1=1;

%% Setting maximal wavelength, change here!
wmax=0.5*n0*(a+b);
ww=wmax; 

while q1==1 

%% Setting initial width of screen
c=10*(ww/a)^2 + (n0+1)*a + (n0-1)*b; 

ww2=[-c/2,c/2];

str=' The wavelength is xx';
str=sprintf(strrep(str,'xx','%g'),ww); 
tt2=text0([.15 .85 .35 .05],str);

q2=1;

while q2==1

%% Display axis
xy=[ww2(1),ww2(2), 0, 1.2];

y0=linspace(ww2(1),ww2(2),N);
fslit='(sign(y1+a/2)+1).*(-sign(y1-a/2)+1)/4'; % shape of one slit
ff=zeros(size(y0));

% Make an equispaced array of slits
for n=1:n0
y2=sprintf('y0+(2*n - n0 - 1)*b/2');
f1=strrep(fslit,'y1',y2);% replace string by other string
f1=eval(f1);
ff=ff+f1;
end

ff=.05*ff;

grill=(2*[1:n0]-n0-1)*b/2; % Center points of array
yy=y0/ww; amp=zeros(size(y0));

for m=1:n0
amp1=fresn(yy-(grill(m)+a/2)/ww)-fresn(yy-(grill(m)-a/2)/ww);
amp=amp+amp1;
end

int=conj(amp).*amp;norm=max(int);int=int/norm;

%% Plot pattern together with array of slits
subplot(1,1,1),
plot(y0,int), axis(xy),hold on,bar(y0,ff), hold off,
string=sprintf('Intensity interference pattern, wavelength = %g',ww);
xlabel('Distance from center of screen. Do you want to zoom in?');
title(string),

rbutt([0.45  0.01  0.15  0.05],'ZOOM IN','q2=1; uiresume');
rbutt([0.6   0.01  0.15  0.05],'PANORAMA','q2=2; uiresume');
gbutt([0.75  0.01  0.15  0.05],'PASS','q2=0; uiresume'); 
uiwait;

if q2==1

%% Zooming in on screen
tt2=text0([.15 .8 .5 .05],' Choose an interval using the mouse!'); 
[ww2,wwy] = ginput(2);
delete(tt2), q1=1;

elseif q2==2
%% Returning to default screen
ww2=[-c/2,c/2];
q1=1;

elseif q2==0 

%% Slider for going to shorter wavelength
tt2=text0([.15 .8 .5 .05],'Change the wavelength usig the slider!'); 
xlabel('Distance from center of screen. Change wavelength with slider!');
obutt,
gbutt([0.75  0.01  0.15  0.05],'PASS','q1=0;uiresume'); 

x=sqrt(ww);

g0 = slider0([.92  .12  .05  .8],[0,sqrt(wmax),x],'q2=0; q1=1; uiresume;');
uiwait; x=get(g0,'Value')+eps; ww=x^2; delete(tt2), delete(g0);

end %% of q2==?

end %% of q2==1

end %% of q1==1
%
rbutt([.45  .01   .15  .05],'REPEAT','uiresume;q1=2;'); 
bbutt([.6  .01   .15  .05],'BACK','close;q1=3;'); 
bbutt([.75  .01   .15  .05],'QUIT','close;q1=4;'); 
uiwait 

% end %%% of q1=1 or q1=2

if q1==3
scatt3d; return

elseif q1==2
q0=1;
obutt,

else
q0=0;
end

end

close,

disp('> Type <diffrac> to start the program again!');
% G Lindblad 970704
