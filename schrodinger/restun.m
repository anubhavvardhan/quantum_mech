
%> The file <restun.m> calculates the transmission coefficient for the
%> 1D Schrodinger equation using multiplication of transfer matrices
%> which are calculated using Airy functions.
%> The potential is chosen to illustrate resonant tunneling.
%> For a general discussion of the physics of resonant tunneling devices 
%> see Capasso & Datta: Physics Today Feb 1990, p74. 
%> For numerics (using other methods) see for instance: 
%> Mendes & Dominguez-Adame: Am J Phys 6282), 143 (1994) 
%> NOTE: for real resonant tunneling devices the one-electron 
%> approximation is insufficient! However, the narrow resonances are 
%> displayed quite effectively by the present simple model.
%>
%> The potential = two square barriers, height V, and there is a bias B.
%> The particle comes in from x = + infinity with energy E, we start 
%> the integration %> from the transmitted wave in x = 0, energy E+B.
%> Variables in thismodel are V (scalar), E,B vectors of same dimensions
%> where you can choose the upper and lower bound. It is recommended
%> that one of them is chosen as constant, otherwise the display will be
%> a bit confusing.
%> At x = L we adapt the incoming and reflected waves and renormalize
%> to find transmission coefficient. 
%> The reflection coefficient is displayed as a function of bias energy.
%> Notation : Standard Schrodinger equation with mass=hbar=1.
%> Theory: Messiah Vol 1, Chapter 3.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 980412

clear;  close;disp('> Welcome to <restun>');

string={' RESONANCES IN QUANTUM TUNNELING'
''
' This file gives a graphic illustration of resonant tunneling '
' of a particle in a potential consisting of two square barriers'
' and a bias voltage.'
''
' For a general discussion of the physics of resonant tunneling  '
' devices see Capasso & Datta: Phys Today February 1990, p74.'
' The general principle is that quasi-bound states in the quantum  '
' well between the two barriers will allow a transmission coefficient'
' close to 1 when the basis voltage is such that the incident energy '
' is resonant with the quasi-bound state.'
''
' Now press the button to display the form of the potential!...........'};

tt1 = text0([.15  .4 .75 .4], string);

cbutt(.75,.01);
delete(tt1),


%%%%%%%%%%%% SOME FIXED PARAMETERS %%%%%%%%%%%%%%%%%%%%%%
L=4;% Length of device % CHANGE HERE
d=1;% Width of barriers	< L/2 % CHANGE HERE

%%%%%%%%%%%% PLOTTING THE POTENTIAL %%%%%%%%%%%%%%%%%%%%%%%%

x=linspace(-3,L+2,500)'; % Space lattice
pot=(sign(x) - sign(x-L/4) + sign(x-3*L/4) - sign(x-L)); % Barriers
pot1=8*pot + x.*(sign(x)- sign(x-L))+L*(sign(x-L)+1); % Just for display
plot(x,zeros(size(x))),axis([-3  L+2 -.5 36]);
hold on 
l1=plot(x,pot1);
set(l1,'LineWidth',2),
axis('off'),
title('Configuration of resonant tunneling system: barriers and bias.'),
l1=text(4.2,16,'Incident wave');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(4.2,14,'Energy = E');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(-3,12,'Transmitted wave');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(-3,10,'Energy = E + B');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(4.2,6.5,'Bias = B');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(-3,22,'Barrier height = V');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(-3,20,'Barrier width = d');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(-3,18,'Device length = L');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(-0.3,-2,'x = 0');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),
l1=text(3.7,-2,'x = L');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),

l1=text(-3,26,'DEFINITION OF PARAMETERS');
set(l1,'FontName','palatino'),set(l1,'FontSize',18),

hold off
xlabel('Press the button to continue!');

%%%%%%%%%%%%% CHOICE OF SOME PARAMETERS %%%%%%%%%%%%%%%



E0=(2*pi/L)^2/2;

M=200; % Number of lattice points for energy/bias.  Choose to suit your computer.  

W=linspace(0,1,M);
cbutt(.75,.01),

string=str2mat(' The parameters chosen so far are:',...
sprintf(' Total length L = %g, width of barriers = %g', L,d),...
' From the value of L follows that the lowest quasi-bound state can be ',...
sprintf(' estimated to be near E0 = %g if V is large and B small.', E0),...
' Now choose the other parameters of the problem!');  

tt1 = text0([.15  .75 .75 .15], string);
cbutt(.75,.01),delete(tt1),

%%%%%%%%%%%% INPUT: ENERGY PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
Q1=3;
			while Q1 < 4

			if Q1==3
str='25';
txt={'Set the barrier height V at right:   V = '};
t2 = text0([.15  .8 .75 .05], txt);
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume'),
VV=edit1([.7 .8 .2  .05],str);
VV=eval(VV);
delete(t2),
			
			Q1=2;
			
			end
			
			if Q1 > 1 

str='1';
txt={'Set the bias as a positive number B > 0, B = '};
t2 = text0([.15  .8 .75 .05], txt);
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume');
BB=edit1([.7 .8 .2  .05],str);
BB=eval(BB);
delete(t2),

str='[0,3]';
txt={'Set the incident energy E as an interval [E0,E1] = '};
t2 = text0([.15  .8 .75 .05], txt);
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume');
EE=edit1([.7 .8 .2  .05],str);
EE=eval(EE);
delete(t2),

			end

E=EE(1)+W*(EE(2)-EE(1)); % Energy lattice

Y = E + BB;          %  Energy of transmitted electrons

Y1=ones(size(Y));

wave1=sqrt(2*E+eps); %  Wave no of incident electron (from right)

wave2=sqrt(2*Y+eps); %  Wave no of transmitted electrons (going left)

wave=[Y1; i./wave1; -i*wave2; wave2./wave1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mm1=trf2([0,d], [VV,VV+BB/4], Y);
mm2=trf2([d,L-d], [BB/4,3*BB/4], Y);
mm1=mult4(mm1,mm2);
mm2=trf2([L-d,L], [VV+3*BB/4,VV+BB], Y);
mm1=mult4(mm1,mm2);
DD=sum(mm1.*wave); TT=2./DD;   % this is the transmitted amplitude

%%%%%%%%% Plotting the results

plot(E,wave2.*abs(TT).^2./wave1 );

title('Transmission coefficient as a function of the energy of incident electrons');

xlabel(' Press one of the buttons');

rbutt([.45  .01  .15  .05],'ZOOM IN','uiresume;Q1=1;'); 
rbutt([.6  .01  .15 .05],'NEW BIAS','uiresume, obutt, Q1=2;');
gbutt([.75  .01  .15  .05],'PASS','uiresume;Q1=6;'); 
uiwait,

if Q1==6

rbutt([.25 .01  .2  .05],'NEW BARRIER','uiresume, obutt, Q1=3;'); 
rbutt([.45  .01  .15  0.05],'REPEAT','uiresume, Q1=4;'); 
bbutt([.6  .01  .15  0.05],'BACK','close;Q1=5;'); 
bbutt([.75  .01  .15  0.05],'QUIT','close;Q1=6;'); 
uiwait;

end 
				
				if Q1==1
				
				tt1= text0([0.2 0.85 .6 .05],'Select a new interval using the mouse!');
				drawnow;
				[EE,Y]=ginput(2);
				delete(tt1),
				end		
				
		end 
				
				if Q1==4
				restun, return,
				elseif Q1==5
				scatt, return,
				end;

disp('> Type <restun> to do this again!');













