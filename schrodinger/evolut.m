
%> The file <evolut.m> calculates the evolution of a wave packet on 
%> the unit interval using an iterative method for the propagation. 
%> Algorithm = Visscher scheme. System = particle in a box with
%> boundary conditions: u=0  at x=0,1 and an additional potential.
%> Initial wave function = Gaussian multiplied by plane wave.
%> Reference: Visscher: Computers in Physics, nov/dec 1991. 


clear; close; disp('> Welcome to <evolut>');

f1=figure;
txt={' The file <evolut.m> calculates the evolution of a wave packet on ' 
' a finite interval using an iterative method for the propagation. ' 
' ' 
' System = Schroedinger particle in [0,1] with u(0) = u(1) = 0  ' 
' and a potential V(x). ' 
' Initial state = Gaussian multiplied by plane wave. ' 
' Algorithm given by Visscher: Computers in Physics, nov/dec 1991. ' 
' ' 
' There are a number of choices for the  potential. ' };

tt1=text0([.15 .3 .75 .5],txt);
t1=title('Evolution of wave packets');
set(t1,'FontSize',18); axis('off');

cbutt(.75, .01); delete(tt1);

%%%%%%%%%%%%% CHOICE OF POTENTIAL %%%%%%%%%%%%

q3=4;

while q3 < 5

	if q3==4

%%%%%%%%%%%%
p1=subplot(2,3,1);
q=get(p1,'Position');
axis([0 1 0 1]);axis('off');%axis('equal')
l1=line([0,0],[0,1]);
set(l1,'LineWidth',3);
l1=line([0,1],[0,0]);
set(l1,'LineWidth',2);
l1=line([1,1],[0,1]);
set(l1,'LineWidth',2);
gbutt([q(1)+0.05  q(2)+ .25  0.1 .05],'#1','uiresume,q1=1;'); 

%%%%%%%%%%%%	
p2=subplot(2,3,2);
q=get(p2,'Position');

axis([0 1 0 1]);axis('off');%axis('equal')
l1=line([0,0],[0,1]);
set(l1,'LineWidth',3);
l1=line([0,1],[0,.5]);
set(l1,'LineWidth',2);
l1=line([1,1],[.5,1]);
set(l1,'LineWidth',2);
gbutt([q(1)+0.05  q(2)+ .25  0.1 .05],'#2','uiresume,q1=2;'); 

%%%%%%%%%%%%
p3=subplot(2,3,3);
q=get(p3,'Position');
axis([0 1 0 1]);axis('off');%axis('equal')
l1=line([0,0],[0,1]);
set(l1,'LineWidth',3);
l1=line([0,.5],[0,0]);
set(l1,'LineWidth',2);
l1=line([.5,.5],[0,.5]);
set(l1,'LineWidth',2);
l1=line([.5,1],[.5,.5]);
set(l1,'LineWidth',2);
l1=line([1,1],[.5,1]);
set(l1,'LineWidth',2);
gbutt([q(1)+0.05  q(2)+ .25  0.1  .05],'#3','uiresume,q1=3;'); 

%%%%%%%%%%%%
p4=subplot(2,3,4);
q=get(p4,'Position');
axis([0 1 0 1]);axis('off');%axis('equal')
l1=line([0,0],[0,1]);
set(l1,'LineWidth',3);
l1=line([0,.2],[0,0]);
set(l1,'LineWidth',2);
l1=line([.2,.2],[0,.5]);
set(l1,'LineWidth',2);
l1=line([.2,.4],[.5,.5]);
set(l1,'LineWidth',2);
l1=line([.4,.4],[.0,.5]);
set(l1,'LineWidth',2);
l1=line([.4,1],[.0,.0]);
set(l1,'LineWidth',2);
l1=line([1,1],[.0,1]);
set(l1,'LineWidth',2);
gbutt([q(1)+0.05  q(2)+ .25  0.1  .05],'#4','uiresume,q1=4;'); 

%%%%%%%%%%%%
p5=subplot(2,3,5);
q=get(p5,'Position');
xx=linspace(0,1);
yy=4*(xx-0.5).^2;
l1=plot(xx,yy);
axis([0 1 0 1]);axis('off');%axis('equal')
set(l1,'LineWidth',2);
gbutt([q(1)+0.05  q(2)+ .25  0.1  .05],'#5','uiresume,q1=5;'); 

%%%%%%%%%%%%
p6=subplot(2,3,6);
axis([0 1 0 1]);axis('off');%axis('equal')
q=get(p6,'Position');

 
t1=suptitle('Choose between the different potentials!');
bbutt([.75  .01 .15 .05],'BACK','uiresume, q1=6;'); 
% bbutt([.6  .01  0.15 .05],'MAIN MENU','close;q1=7'); 
% bbutt([.75  .01  0.15 .05],'QUIT','close;q1=8;'); 

uiwait;
  
		 if q1==1
		 potential='zeros(size(X))';
y=str2mat('> A free particle in a potential box moves under the kinetic Hamiltonian.',...
		'> We can see the dispersion of the wave packet as it bounces between the walls',...
		'> of the box and eventually how it spreads all over the box and reduces to',... 
		'> random waves moving about. '); 	 

		 elseif q1==2
		 potential='q2*X'; 
y=str2mat('> The potential is linear V(X) = k X. You have to choose the constant k',...
		'> and you had better choose it of the order of 10^4 or more. The initial wave ',...
		'> packet can be chosen to be at rest near the top of the hill. The near harmonic',...
		'> nature of the spectrum means that the wave packet is reassembled in an almost ',...
		'> periodic fashion.');

		 elseif q1==3
		 potential='q2*0.5*(1+sign(X-0.5))';
y=str2mat('> The potential is a potential step in a finite box. ');

		 elseif q1==4
		 potential='q2*0.5*(sign(X-0.2) + sign(0.4-X))';
y=str2mat('> The potential is a barrier or hole in a finite box .'); 

		 elseif q1==5
		 potential=' q2*4*(X-0.5).^2 '; 
y=str2mat('> A harmonic potential in a finite box, centered at X=0.5.',...
		'> If the multiplier is large enough we will get a spectrum close to ',...
		'> harmonic for the lowest level. Choose multiplier of the  order of 10^4.');
		 
		 elseif q1==6
		 wavepac;return;
		 
		 elseif q1==7
		 start; return;

		 elseif q1==8
		 close; return;
		 end
 	
	 	 disp(y);	

 	q3=3;

	end %% of q3==4
		
	if q3==3
	
txt={'Choose a multiplier of'
	'order 10^4 (except for #1),'
	'the size of the space grid,'
	'a starting point '
	'a starting velocity '
	'and the number of pictures'};

tt1=text0([.65 .1 .3 .35], txt);

		if q1==1
		q2=0; 
		else
		
tt2=text0([.65 .15 .3 .05],'Enter the multiplier');
str='10^4';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
ee=edit1([.8 .1 .15  .05],str),
delete(tt2),
q2=eval(ee);

		end

tt2=text0([.65 .15 .3 .05],'Enter the # grid points');
str='500';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
ee=edit1([.8 .1 .15  .05],str); 
delete(tt1,tt2),
N=eval(ee);

%%% SHOW POTENTIAL %%%

	X=linspace(0,1,N); %
	V=eval(potential);
	maxV=max(V);
	minV=abs(min(V))+eps;

%%% CONSTRUCTION OF THE HAMILTONIAN %%%%%%%%

uni=speye(N);

L = 1;              % system extends from 0 to L 
h = L/(N-1);        % grid cell size
hh=-1/(2*h^2); 	% factor in front of Laplacian

taumax=min(2/minV,2/(maxV + 2/h/h));
disp(sprintf('> The maximal time step is:%g ',taumax));

VV=sparse(1:N,1:N,V);
od=ones(N-1,1); dd=ones(N,1); 
DD=[sparse(1:N-1,2:N,1);sparse(1,1:N,0)]+[sparse(2:N,1:N-1,1), sparse(1:N,1,0)]-2*speye(N);

%tau = input('> Enter the time step  >>>  ');
tau=0.9*taumax;  HH=tau*(hh*DD+VV);

%%%%%%%%%%%%%

	q3=2; 

	end %% of q3==3

	
	if q3==2

%%% INPUT OF WAVE PACKET PARAMETERS %%%
txt={'Enter the center of the wave'
	'packet, a number in (0,1)'};
tt2=text0([.65 .1 .3 .3],txt);
str='.8';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
ee=edit1([.8 .1 .15  .05],str); 
delete(tt2),
x0=eval(ee);

txt= sprintf('Choose a velocity smaller than %g ',N/L);
tt2=text0([.65 .1 .3 .3],txt);
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
str='20';
ee=edit1([.8  .1 .15  .05],str);
k0=-eval(ee); delete(tt2),

sigma0 = L/30;   % Standard deviation of the wavefunction
Norm = sqrt(h)/(sqrt(sigma0*sqrt(pi)));  % Normalization

txt= {'Enter the number of pictures'};
tt2=text0([.65 .1 .3 .3],txt);
str='100';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
ee=edit1([.8  .1 .15  .05],str);
M=eval(ee);

Y = Norm * exp(i*k0*X').*exp(-(X'-x0).^2/(2*sigma0^2));
Y(1) = 0;  Y(N)=0;

Z0 =real(Y.*conj(Y));   % Record initial condition
NZ0=sum(Z0); close,
plot(X,Z0/NZ0);xy=axis;
title('Initial wave packet (probability density)');
xlabel('Press button to continue!');
cbutt(.75,.01);
rbutt([.75 .01 .15 .05],'WAIT..',''), drawnow;

%%% INITIAL VALUES
R0=real(Y);
% Calculate second initial value (using Crank-Nicholson): 
WW=uni-i*HH*0.125;
Y=WW*Y; Y=WW'\Y; Y=WW*Y; Y=WW'\Y;

I0=imag(Y); 
clear WW DD 

%%% STATISTICS OF STATE
average=[]; deviation=[];
RR=R0; II=I0;

%%% INPUT OF GRAPHICS PARAMETERS %%%
kmean=sqrt(maxV + k0^2);

Q=floor(1/40/kmean/tau);
disp(sprintf('> Have chosen %g iterations per picture',Q));
disp(sprintf('> Now performing  %g iterations',M*Q));

M0=1;	q3=1;
	
	end
	
	if q3==1

	for  m=M0:M

		for  q=1:Q

		RR=RR+HH*II;		II=II - HH*RR;

		end

	Z1=RR.^2 + II.*(II+HH*RR);
	plot(X,Z1);	axis(xy); axis('off');
	title(sprintf(' Picture number %g ',m));
	drawnow;
	
	average=[average,X*Z1];
	deviation=[deviation,sqrt(X.^2*Z1-(X*Z1)^2)];

	end

cbutt(.75,.01);

M0=M+1;MM=1:M; 
plot(MM,[average;deviation]);
title('Average <x> and standard deviation as functions of time');
% 

		end %% of q3==1


q3=6;


rbutt([.55  .01  .2  .05],'More pictures','uiresume, q3=1;');
cbutt(.75, .01);


if q3==6

rbutt([.55  .01  .2  .05],'New velocity','uiresume, q3=2;');
cbutt(.75, .01);

end 

if q3==6

rbutt([.55  .01  .2  .05],'New multiplier','uiresume, q3=3;'); 
cbutt(.75, .01);

end 

if q3==6

rbutt([.4  .01  .2  .05],'New potential','close; q3=4;');
bbutt([.6  .01 .15  .05],'BACK','close;q3=5;'); 
bbutt([.75  .01   .15  .05],'QUIT','close;q3=6;');
uiwait;

end 



		if q3==1
		
txt=sprintf('Choose a number of pictures > %g ',M);
tt2=text0([.3 .2 .6 .05], txt);
str='200';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
ee=edit1([.75 .2 .15 .05],str);delete(tt2),
M1=eval(ee);


		M0=M+1;M=M1;
		
		elseif q3==5
		wavepac;return;
		end

end %%% of q3

disp('> Type <evolut> to do this again. Bye!')


