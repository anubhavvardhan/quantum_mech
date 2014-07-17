
%> The file <trix.m> calculates the energy as a function of Bloch wave number
%> for a periodic 1D discretized  Schrodinger equation.
%> It also displays the allowed and forbidden energy bands.
%> The Hamiltonian has a diagonal part which is a random vector of variable 
%> strength and constant next-to-diagonal elements allowing 'hopping' between
%> sites (this amplitude is normalized to 1).
%> The whole forms an infinite Jacobi matrix.
%> There is a period which can be chosen, inside the unit cell the diagonal 
%> elements are independent random variables with a rectangular distribution.
%> The transfer matrix of the unit cell is calculated for a range of energies 
%> and hence the Bloch momentum as a function of energy.
%> CHOICE: period (integer) and  relative strength of the diagonal part.
%> See how the allowed bands go to zero with increasing strength, the number 
%> of allowed bands is equal to the period.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125

clear;close; disp('> Welcome to <trix>!');

txt={' DISCRETE PERIODIC QUANTUM MODEL' 
' ' 
' This program solves a discrete quantum model of a periodic system.' 
' The configuration space is Z (the integers). The Hamiltonian is an' 
' infinite Jacobi matrix. The diagonal is periodic, with a choice of' 
' period and a multiplier. The values in the unit cell are chosen ' 
' using the random number generators in MATLAB.' 
' The next-to-diagonal elements are all equal to 1. ' 
' ' 
' This program calculates the energies of the Bloch solutions as' 
' functions of Bloch wave number and the allowed energy bands.' 
' The method is based on the calculation of the transfer matrix.'};

tt1=text0([.15 .5 .75 .35], txt);

subplot(2,2,3);
N=16;
z0=[0,1,zeros(1,N-2)];
m0=toeplitz(z0);
spy(m0,30,'b'), %axis('off'),
hold on,
z0=[1,0,0,0];
z1=[z0,z0,z0,z0];
spy(diag(z1),30,'r'), hold on,
z0=1-z0;
z1=[z0,z0,z0,z0];
spy(diag(z1),30,'g'), hold off,
cbutt(.75, .01), delete(tt1),

q0=1;

while  q0==1|q0==2

if q0==1

tt1=text0([.15 .5 .75 .1], 'Choose length of unit cell, integer > 1: ');
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='5';
ee=edit1([.75 .5 .15 .05],str);
q1=eval(ee); delete(tt1), 

vv=rand(1,q1); %random vector for potential
avv=sum(vv)/q1;
vv0=2*(vv-avv); %random potential with factor q2, chosen to have average 0
norm=sum(vv0.^2);
vv0=vv0/sqrt(norm);

subplot(2,2,4),
bar(vv0),
title('Random potential (normalized) '),
xlabel('Press the button to continue!');
cbutt(.75,.01);

q2=1;
disp('> Default multiplier = 1');
string0='Multiplier = q2';
end %%%%% of q0=1

%% now q0=2

string=sprintf(strrep(string0,'q2','%g'),q2);  

vv1=q2*vv0; %random potential with factor q2, chosen to have average 0

%%%%%%%%%%%%% Choice of energy range 
emax=1.2+ .5*q2;
ee=emax*linspace(-1,1,500); %% Change here!!
e1=ones(size(ee));
e0=zeros(size(ee));
m0=diag([1,1]);
%%%%%%%%%%%%% Replace 2x2 matrices by 4-vectors 
%M0=[1,1,0,0]';
M0=[1;0;0;1]; %%% [M11; M12; M13 ;M22]
M=M0*e1;

%% Repeated matrix multiplication:

		for n=1:q1

MI=[2*(vv1(n) - ee); -e1; e1;e0];
M=mult4(MI,M);

		end
		
MM=0.5*(M(1,:)+M(4,:)); % 0.5*(trace of the transfer matrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,1,1),
MM1=abs(MM);
MM1=crop(MM1,10000000, 1);
MM1=acosh(MM1+eps); %MM1=crop(MM1,5,0);
plot(ee,MM1), axis([-emax emax 0 5]),legend(string);
title('Positive Liapounov exponent as a function of energy');
xlabel('The exponent is zero in allowed, positive in forbidden region');
cbutt(.75,.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MM0=crop(MM,10,-10);
plot(ee,MM), axis([-emax emax -5 5]), hold on,
plot(ee, [e1;-e1],'.- r'),hold off,
title('Trace of transfer matrix as a function of energy'),
xlabel('Allowed values in range [-1,1],'),
cbutt(.75,.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MM=crop(MM,1,-1); %% Restricts the values of MM to [-1,1].
plot(ee,MM),axis([-emax emax -1 1]),
title('Normalized trace of transfer matrix  for allowed energy values'),
xlabel('The Bloch wave number is the inverse cosine of this function'),
y=acos(MM); % Bloch wave number.
w2=[-0.1,pi,min(ee),max(ee)];
MM0 = 1-abs(MM); % This is zero where MM has values +/- 1.
MM0=~MM0; % = 1 where MM has values +/- 1, else zero.
cbutt(.75,.01),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(y,ee);axis(w2);
xlabel('Bloch wave number K*L. Change multiplier using the slider!');
ylabel('Energy of Bloch solution and allowed bands');
title('Energy of Bloch waves versus wave number, forbidden bands marked');
hold on; barh(ee,-0.1*MM0); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q0=2;
gbutt([.75  .01  .15 .05],'Pass','uiresume, q0=1;'); 
q2 = slider1([.92  .12  .05  .8],[0,2,q2]);

if q0==1


rbutt([.3  .01  .15   .05],'New V(X)','uiresume; q0=1;');
bbutt([.45  .01  .15  .05],'BACK','close;q0=3;'); 
bbutt([.6   .01  .15  .05],'MAIN MENU','close;q0=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q0=0;'); 
uiwait;

if q0==1
obutt,
end

end

	end %%%% of q0 iteration

	if q0==3
	per1d;	return;
	elseif q0==4
	start; return;
	end 

close;
disp('> Type <trix>  to do this again ');

