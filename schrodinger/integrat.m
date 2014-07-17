
%> The file <integrat.m> displays some basics on the numerical integration
%> of the 1D Schroedinger equation.  First the lack of numerical stability
%> for integration in the classical forbidden region is displayed and the 
%> differences in stability of the different algorithms indicated.
%> The effects of the instability on eigenvalues and eigenfunctions is 
%> considered in a simple example.
%> © Goran Lindblad - gli@theophys.kth.se

clear; close; format long;
global potential energy
disp('> Welcome to <integrat>');

q1=1;

while q1==1 

t1= title('Integration of the Schroedinger equation');
set(t1,'FontSize',18);
axis('off');


txt= {' '
' The program <integrat.m> displays some basics on the numerical' 
' integration of the 1D Schroedinger equation. We write the DE as'
' '
'    D^2 u(x) = 2 [ V(x) - E ] u(x),   a <= x <= b'
''
' Starting from  x = a, we need two initial values to perform the '
' integration, for instance u(a) and Du(a). If the potential is '
' sufficiently regular there will then be a unique solution in [a,b].'
' Each such solution is a linear combination of two fundamental '
' solutions:'
''
'    u(x):  u(a) = 1,  Du(a) = 0,'
''
'    v(x):  v(a) = 0,  Dv(a) = 1.'
''
' We can calculate these solutions for several choices of the '
' energy  E. We find that by variation of E it is possible to fulfill a '
' boundary condition in x = b. Here we will consider BCs of the form  '
'  A u(b) + B v(b) = 0. In particular we will calculate v(x) for '
' different values of E and show how we can find the discrete set of '
' energy values for which v(b) = 0. The integration is performed using'
' a standard MATLAB routine, or using the Numerov method.'
' ' 
' Now press the button to continue!' };

t1=text0([.12 .15 .78 .7],txt);
cbutt(.75, .01), delete(t1),

energy=21; len=1; vmax=100; emin=15; emax=2*vmax;  q=1;
xx=[0,len]; y0=[0,1];

potential='Y*sin(pi*(X - .5)).^2';

potential=sprintf(strrep(potential,'Y','%g'),vmax);

X=linspace(0,len); Y=eval(potential);
l1=plot(X,[Y;ones(size(X))*energy]);
title('This is the potential V and the chosen energy value!');
set(l1,'LineWidth',2);

[X1,Y1]=ode23s('sde',xx,[1,0]);
[X2,Y2]=ode23s('sde',xx,[0,1]);

cbutt(.75, .01);
plot(X1,Y1); hold on,
plot(X2,Y2,'+'); hold off,
title('The fundamental solutions u (--), v (++), with their derivatives Du, Dv. ');
xy=[0 len -100 +100];
axis(xy);

txt=' This is v(X) as a function of energy when v(0) = 0, Dv(0) = 1.';
txt=sprintf(strrep(txt,'X','%g'),len);
ee=linspace(emin,emax); n2=1000;
ww=numerov1(0,len,n2,potential,1,ee,y0(1),y0(2)); 

cbutt(.75, .01);
plot(ee,ww(1,:)), hold on;
title(txt); xlabel('Energy');

w0 = findzero(ee,ww(1,:));
disp(w0');

cbutt(.75, .01);
plot(w0,zeros(size(w0)),'or'); hold off;
title('The energies for which v = 0 at both endpoints of the interval'); 
cbutt(.75, .01);

while q==1

[X1,Y1]=ode23s('sde',xx,y0);
nn=length(X1); n2=5*nn;

txt='The energy = E, change with slider until v = 0 in right end point!';
txt=sprintf(strrep(txt,'E','%g'),energy);
xy=[0 len -1 1];
plot(X1,[Y1(:,1), .2*Y1(:,2)]);%axis(xy);
tt1= text0([.15  .15 .75 .03], txt);
title('The solution v(X) and 0.2 x Dv(X) for variable energy');
xlabel('Use Delete to delete all except last data!');
Z1=[Y1(nn,:), Y1(nn,2)/Y1(nn,1)];

add=0;
g0 = slider0([.91  .15  .05  .75],[-1,1,add],'q=1;hold on; uiresume;');

rbutt([.45  .01  .15  .05],'Delete','hold off, q=1;uiresume;');
rbutt([.6   .01  .15  .05],'Numerov','hold on, q=2;uiresume;');
gbutt([.75  .01  .15  .05],'CONTINUE','hold off, q=3;uiresume;');
uiwait; 

if q==1

add = get(g0,'Value'); energy =energy*(1 + .5*sign(add)*add^2);
delete(g0);

elseif  q==2

obutt, delete(g0,tt1);
Y2=numerov2(xx(1),xx(2),n2,potential,1,energy,y0(1),y0(2)); 
X2=linspace(xx(1), xx(2),n2);
plot(X2,Y2','r'); hold off;
txt='Note the possible divergence of the numerical solutions!!';
tt2=text0([.15 .12 .75 .03], txt);;
gbutt([.75 .01  .15  .05],'CONTINUE','hold on, q=1;uiresume;');
uiwait; delete(tt2);

end

end

close

if q==3

txt= {' '
' THE NUMERICAL INSTABILITY'
''
' The numerical integration of the Schroedinger equation meets with'
' the problem of numerical instability. In a classically forbidden'
' region where E < V(x) the two linearly independent solutions will'
' both be exponentially increasing in at least one direction.' 
' The effects of this instability on eigenvalues and eigenfunctions '
' will be shown in a simple example.' 
' ' 
' We use a potential V(x) in [0,5] which is a positive constant V_0 '
' in  [0,3] and zero in [3, 5]. Choose a value for the energy E such '
' that E - V_0 < 0, hence making [0,3] into a forbidden region.' 
' ' 
' We start the integration from x = 0 and choose values of {u, Du} '
' there which single out the exponentially damped solution. Due to ' 
' the numerical instability the exponentially diverging solution  ' 
' will eventually creep in. The integration is performed using the ' 
' algorithms: ode45, ode23 and numerov2. ' 
' ' 
' Now press the button to continue!' };


t1=text0([.12 .15 .78 .7],txt);
% title('The problem of numerical instability');

%% We define a potential, which is positive on an interval [x0,x1], then
%% zero in [x1,x2].
P = ' V*0.5*(1 + sign(x1 - X))'; 

%% Some parameters: %% You can change here
x0=0; x1=3; x2=5; V=100; 

P=sprintf(strrep(P,'x1','%g'),x1); %% Replace string by value 

txt='> The chosen potential is = '; 

disp([txt,P]);

txt='> The value of V = VV ';  y=sprintf(strrep(txt,'VV','%g'),V); 

disp(txt);

cbutt(.75,.01), delete(t1),

X=linspace(x0,x2);
Y=eval(P);
l1=plot(X,Y); axis([x0 x2 -0.2*V, 1.2*V]); 
set(l1,'LineWidth',3); 
title('This is the given potential');
xlabel('Press the button to continue!');

txt={' We will integrate the Schroedinger'
	' equation from 0 to 5 at an energy '
	' 0 < E < V_0 = 100, consequently the  '
	' interval [0,3] is classically forbidden.' };
	
t1=text0([.15 .25 .4 .15], txt);
cbutt(.75,.01);
delete(t1);
 
txt={'Choose a value  0 < E < Vmax using '
	'the slider. The instability is more '
	'obvious the lower the energy!'};
t1=text0([.15 .25 .4 .15],txt);
	
q=slider1([.92  .12  .05  .8], [0,1,.5]);
	
rbutt([0.75  0.01  0.15  0.05],'WAIT..','');
drawnow, delete(t1); % delete(g0);
energy=q*V;

if isempty(energy)
energy=V/2;
end

potential=sprintf(strrep(P,'V','%g'),V); %% Replace string by value

X=0; V=eval(potential);

kappa=sqrt(2*(V-energy));	%% Wave number
Y0=1;	DY0=-kappa*Y0;		%% Starting conditions in x0
y0=[Y0;DY0];
xy=[x0, x2, -1.5, 1.5];
xxx=[x0,x2];

[X1,Y1]=ode45('sde',xxx,y0);
	
plot(X1,Y1(:,1),'+r'),
title('Decaying solution of the differential equation'), axis(xy), 
drawnow, hold on;
 
[X2,Y2]=ode23('sde',xxx,y0);

plot(X2,Y2(:,1),'og'); 
drawnow,

N=500; %% The number of integration points in Numerov scheme

X=linspace(x0,x2,N); 

ww=numerov2(x0,x2,N,potential,1,energy,Y0,DY0); 

plot(X,ww);
txt={'The Numerov scheme is the least stable!'};
t2=text0([.15 .21 .5 .05],txt);
txt={'Solutions: ode45=+++, ode23=ooo, numerov2 = ---'};
t1=text0([.15 .15 .6 .05],txt);
hold off,

cbutt(.75,.01);
delete(t2),

txt={'We now do the same thing but single out the exponentially  '
 ' increasing solution, identical to the diverging one.  '
 ' The integration is performed using the same algorithms:'
 ' ode45, ode23 and numerov2'};

t2=text0([.15  .21  .65 .12],txt);
cbutt(.75,.01), delete(t2),
rbutt([0.75  0.01  0.15  0.05],'WAIT..',''),
drawnow;


Y0=1;	DY0=kappa*Y0;		%% Starting conditions in x0
y0=[Y0;DY0]; xxx=[x0,x2];

[X1,Y1]=ode45('sde',xxx,y0);

plot(X1,Y1(:,1),'+r');
title('The divergent solution'), drawnow, hold on;
 
[X2,Y2]=ode23('sde',xxx,y0); 
plot(X2,Y2(:,1),'og'); drawnow

N=500; %% The number of integration points in Numerov scheme

X=linspace(x0,x2,N); 

ww=numerov2(x0,x2,N,potential,1,energy,Y0,DY0); 

plot(X,ww);

txt={'The Numerov scheme now works as well as the others!'};

t2=text0([.15 .21 .65 .05], txt);

hold off, cbutt(.75,.01), delete(t1,t2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

txt={' '
' FINDING THE EIGENVALUES '
''
' We now return to the problem of finding energy levels by numerical'
' integration of the Schroedinger equation in an interval [0,1].'
' The potential V is first chosen to be constant, so we know the  '
' exact solutions to be harmonic functions. '
'' 
' Choose the BCs  v(0) = v(1) = 0.'
''
' We use the Numerov algorithm to integrate from x = 0 to x =1 for a '  
' finite set of energy values, of the order of 100 - 500 points. '
' We only need the solution in x = 1 as a function of E.'
''
' The values of E for which v(x=1,E) = 0 is then estimated using '
' spline interpolation, '
''
' Now press the button to continue!'};


close

t1=text0([.12 .15 .78 .7],txt);
rbutt([0.75  0.01  0.15  0.05],'WAIT..',''), drawnow;

x0=0;	x1=1; N=300;  

q1=[1,0];

X=linspace(x0,x1,N); X0=zeros(size(X));

E=linspace(0,200,300); E1=ones(size(E));

P='ones(size(X))';

DY0=E1; Y0 = q1(2)*E1;

W=numerov1(x0,x1,N,P,0,E,Y0,DY0);

ww=q1(1)*W(1,:)+ q1(2)*W(2,:); %% BC: ww  = 0;
cbutt(.75,.01); delete(t1), 

plot(E,[ww;X0]); xy=axis;
txt={' This is v(x=1,E), the boundary value in x = 1 as a function '
' of the energy E, calculated for a discrete set of values of E.'
' The energy eigenvalues are given by the crossings with the '
' y = 0 axis. Their values are estimated by the algorithm  '
' <findzero> which uses spline interpolation! Then the '
' corresponding eigensolutions are found by integration. '
' '};
t1=text0([.2 .65 .65 .25],txt);
title('Solving the eigenvalues from the boundary conditions'),
ylabel('v(x=1, E)');
xlabel('The  energy parameter E'),
drawnow, hold on;

zz=findzero(E,ww); disp(zz');
cbutt(.75,.01);

plot(zz,zeros(size(zz)),'or'); hold off;

Y0=ones(size(zz));

E=zz; E1=ones(size(E));

if q1(2)==0

DY0=E1; Y0 = q1(2)*E1;

else

Y0=E1; DY0= -q1(1)*E1/q1(2);

end

W=numerov2(x0,x1,N,P,0,zz,Y0,DY0);

cbutt(.75,.01);
delete(t1);


plot(X,W); title('The eigensolutions (non-normalized)');
xlabel('Press the button to continue!');
cbutt(.75,.01); close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

txt={' '
' NUMERICAL INSTABILITY ..'
''
' We now look at the same problem when there is a forbidden interval'
' in the range of integration.  The potential is that of the harmonic'
' oscillator, so we know the exact solutions. '
''
' We again use the Numerov algorithm to find the solutions at the  '
' right endpoint, and from these data estimating the eigenvalues.'
' The eigenfunctions are again integrated using the Numerov algorithm. '
'' 
' We will see that the eigenfunctions are more sensitive than the'
' eigenvalues to the instability. '
''
' Now press the button to continue!'};



t1=text0([.12 .15 .78 .7],txt);
K=32; %%  the constant in harmonic potential

P='K*X.^2'; %% Harmonic potential 

P=sprintf(strrep(P,'K','%g'),K);

x0=-3; x1=3; %% End points 

N=100; %% The number of integration points per unit interval

X=linspace(-3,3,6*N); Y=eval(P); 

cbutt(.75,.01);

delete(t1),
plot(X,Y); 



txt ={' Finding the eigenvalues using <findzero.m>.....'};
% t1=text0([.25 .85 .5 .05],txt);
rbutt([0.75  0.01  0.15  0.05],'WAIT..','');
drawnow;	hold on;

E=linspace(0,150,300); 

Y0=ones(size(E)); DY0=-2*x0*Y0;

W=numerov1(x0,x1,6*N,P,1,E,Y0,DY0);

wronski=W(1,:).*(2*x1*Y0) - W(2,:); %% DY has opposite sign

zz=findzero(E,wronski); %% Finding the eigenvalues

disp(zz'); %% Displaying the eigenvalues

zz1=zz'*ones(size(X));

plot(X,zz1,'r');

title('Energy levels superimposed on potential');

txt={' You can check that the HO spectrum is well approximated!'};
t1=text0([.2 .82 .65 .05], txt);

drawnow; hold off;
cbutt(.75,.01);
delete(t1);
rbutt([0.75  0.01  0.15  0.05],'WAIT..',''), drawnow;
n0=length(zz);

nn=[1:n0];

omega=sqrt(2*K); %% frequency of HO

plot(nn,zz,'or');
axis([1,n0,0,max(zz)+1]);
title('Energy versus order of eigevalue'); hold on 
plot(nn,omega*(nn-0.5),'+'); hold off
xlabel('Numerical values = o, exact harmonic spectrum = +');
txt={' Now finding the eigenfunctions..'};
t1=text0([.2 .8  .4 .05],txt);
rbutt([0.75  0.01  0.15  0.05],'WAIT..','');

Y0=ones(size(zz)); DY0=-2*x0*Y0;

W2=numerov2(x0,x1,6*N,P,1,zz,Y0,DY0);

WW2=W2'*W2;	WW3=diag(WW2);	WW3=1./sqrt(WW3);

WW3=diag(WW3);	W3=WW3*WW2*WW3;	W20=W2*WW3;

delete(t1);

cbutt(.75,.01);
p1=plot(X,W20); 
txt={' The integration fails to get the eigenfunctions because the'
 ' exponentially divergent solution dominates everything!'};
t1=text0([.2 .7  .65 .1],txt);
xlabel('The space coordinate');

cbutt(.75,.01);
delete(t1);


txt={' '
' The instability neutralized....'
''
' We now do the integration in such a way that the relevant '
' solution is always the divergent one.  Start from the left hand '
' end of the interval and integrate to the midpoint x=0. Use that  '
' the solutions are either even or odd to get the BC at x=0.'
' For non-symmetric potentials we start another integration  '
' from the right end point, integrate to x=0, and identify the '
' logarithmic derivatives in this point.'
''
' Integrating using the Numerov algorithm...'};


t1=text0([.15 .15 .75 .7],txt); drawnow;
Y0=ones(size(E));	DY0=-2*x0*Y0;

W=numerov1(x0,0,3*N,P,1,E,Y0,DY0);

ww=W(1,:).*W(2,:); %% We get all even and odd solutions setting ww=0!

zz=findzero(E,ww); % disp(zz'); %% The eigevalues are approx the same!!

Y0=ones(size(zz)); DY0=-2*x0*Y0;

W2=numerov2(x0,0,3*N,P,1,zz,Y0,DY0);

WW2=W2'*W2; WW3=diag(WW2); WW3=1./sqrt(WW3);

WW3=diag(WW3); W3=WW3*WW2*WW3; W20=W2*WW3;

X=linspace(-3,0,3*N); 

cbutt(.75,.01),delete(t1), plot(X,W20);
title('QHO eigenfunctions, displayed on half interval'),
txt={' We now get reasonable approximations of the eigenfunctions!'
     ' The accuracy of the eigenvalues is approximately the same!'};
t1=text0([0.15 .85 .75 .06],txt);


rbutt([.45  .01  .15 .05],'REPEAT','close; q1=1;'), 
bbutt([.6  .01   .15  .05],'MAIN MENU','close;q1=2;'), 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=3;'), 
uiwait; 
end
	if q1==2
	
	close; start; return;
	
	elseif q1==3
	
	q1=0;
	
	end

	end

disp('> Type <integrat> to do this again.  Bye!');


%%% © Goran Lindblad 1996
