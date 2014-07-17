%> <well2.m> solves for the eigenvalues of a quadratic and a circular 
%> well in 2D, of the same area and depth, and compares the spectra
%> Variable parameters, A = area, D = depth.
%> The search for eigenvalues i automatic, using <findzero.m>
%>
%> © Goran Lindblad - gli@theophys.kth.se

%% GL 990428

clear, close; disp('> Welcome to <well2>');

txt={' BOUND STATES OF A 2D WELL ' 
' ' 
' The script <well2.m> solves for approximate eigenvalues for a 2D ' 
' quantum well of quadratic or circular shape of given area A and  '
' constant depth D.' 
''
' It solves for the eigenvalues from the boundary condition, where'
' the solutions regular inside the boundary are continuously joined'
' to the exponentially decaying solutions outside. '
''
' The quadratic case is solved by first solving the 1D problem, and '
' then adding the resulting spectra to the spectrum of the 2D case.'
' In the circular well the radial functions are Bessel functions of'
' types J_m inside, and K_m outside the boundary. The numerical '
' search for eigenvalues proceeds for increasing m = 0,1,... until'
' no more eigenvalues are found. '
''
' You must choose the values of A and D. If you choose the values'
' too large the search for the eigenvalues may fail for the circular'
' case!' 
' '
' Now insert new values of [A,D] in the box >> '
' or press the continue button!'};

%% Constants, CHANGE HERE
nn=1000; % lattice points for energy scale

tt1=text0([.12 .3 .78 .6], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='[ 30 , 5 ]';
ee=edit1([.75 .3 .15  .05],str);
AD=eval(ee); 
A=AD(1); D=AD(2);
rdim=sqrt(A/pi);xdim= sqrt(A);

delete(tt1);

%% square well 
%% xx=linspace(-xdim/2, xdim/2);
%% wave number inside
kk=linspace(eps, sqrt(2*D)-eps,nn);
%% exponential decay outside
kk0= sqrt(2*D - kk.^2);
kk0=real(kk0);

%% even solutions  cos(kk*X) and exp(-kk0*X)
wron= real(- kk.*sin(kk*xdim/2) + kk0.* cos(kk*xdim/2));
plot(kk,wron); axis tight,
title('Finding the eigenvalues for a 1D square well, even and odd solutions'),
hold on,
z1=findzero(kk,wron);
plot(z1,zeros(length(z1)), 'o r')

%% odd solutions  sin(kk*X) and exp(-kk0*X)
wron= real(kk.*cos(kk*xdim/2) + kk0.* sin(kk*xdim/2));
plot(kk,wron,'k');
z2=findzero(kk,wron);

%% there is a false solution kk=0; remove it:
if z2(1)==0

z2=z2(2:length(z2));

end

plot(z2,zeros(length(z2)), '^ r')
hold off

zz=[z1,z2];
zz=sort(zz);
z2=zz.^2/2;
ee= z2 - D;
%

cbutt(.75,.01);

xy=[.5 , length(ee)+.5 , - D , 0];
bar(ee); axis(xy),
%length(ee);
title('Bound state energies for the 1D well');

cbutt(.75,.01);

[e1,e2]=meshgrid(z2,z2);
eee=e1+e2-D;
%% we have to find the negative values in eee!
eee=reshape(eee,length(z2)^2,1);
eee=sort(eee);
eee=crop(eee,0,-D); %% positive eigenvalues set to zero
ind=find(eee);
eee=eee(ind); %% finds nonzero (strictly negative) ones
xy=[.5 ,length(eee)+.5, - D, 0];
bar(eee); axis(xy)

title('Bound state energies for the square 2D well');
length(eee);
cbutt(.75,.01)
ll=length(eee);
[x,y]=stairs(eee);
x=[0;0;x-1;ll;ll];y=[-D;-D;y;eee(ll);eee(ll)];
p1=plot(y,x); axis([-D, 0, 0, ll]),
title('Distribution function for eigenvalues, quadratic case'),
set(p1,'LineWidth',3),
cbutt(.75,.01);
rbutt([.75  .01  .15  .05],'WAIT..','');
drawnow;


%% Now the circular well
%% the energy for magn quantum number m is kmn^2/2 - D
%% kmn is the solution of the boundary condition at r = rdim
%% solution inside is J_m(kk r), outside K_m(kk0 r)
%% D J_m /J_m = m  /z - J_m+1 /J_m
%% DK_m /K_m = m /z - K_m+1/K_m

iter = 1; m=0; ee2=[];

while iter ==1

wron = kk.*besselj(m+1,kk*rdim).*besselk(m,kk0*rdim)./(kk0+eps)...
- besselk(m+1,kk0*rdim).*besselj(m,kk*rdim);

wron=real(wron);
wron=crop(wron,1.e+10, -1.e+10);

 
z1=findzero(kk,wron);

 
 if isempty(z1)
 
 iter = 0;
 
 else
 
 if z1(1) < eps
 
 z1=z1(2:length(z1));
 
 end
  
z2=z1.^2/2;  ee1= z2 - D; 
ee1=crop(ee1,0,- D);
ind=find(ee1);
ee1=ee1(ind);

if m > 0

ee1=[ee1,ee1]; % double degeneracy!!!
ee1=sort(ee1);

end
 
ee2=mad(ee2,ee1);
m=m+1;

if isempty(ee1)
iter = 0
end 
 
 end 
 
 end 
 
emin= - max(max(-ee2))
le=size(ee2);
bar([0:le(1)-1],ee2,'b');
axis([-.5 , le(1)-.5, emin, 0]);
title('Circular case: bound state energy spectra for  m = 0, 1, ..');
xlabel('Angular momentum m');
cbutt(.75, .01)

ee2=reshape(ee2,le(1)*le(2),1);
ee2=sort(ee2);
bar(ee2);axis(xy),
title('Bound state energies of circular well');
l2=length(ee2);

cbutt(.75,.01)


[x2,y2]=stairs(ee2);
x2=[0;0;x2-1;l2;l2;l2;l2];y2=[-D;-D;y2;ee2(l2);ee2(l2);0;0];
p1=plot(y2,x2); axis([-D, 0, 0, ll]),
title('Distribution function for eigenvalues, circular case'),
set(p1,'LineWidth',3),
cbutt(.75,.01);
hold on,
p1=plot(y,x,'r');
set(p1,'LineWidth',3),
title('Distribution function for eigenvalues, circular and quadratic cases'),
hold off,

rbutt([.3  .01  .15 .05],'REPEAT','close; q1=1;'),
bbutt([.45  .01   .15  .05],'BACK','close;q1=2;') 
bbutt([.6  .01   .15  .05],'MAIN MENU','close;q1=3;'), 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=4;'), 
uiwait; 

if q1==1
well2;
elseif q1==2
close; boundst; return;
elseif q1==3	
close; start; return;
else
close;
disp('> Type <well2> to do this again!');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
