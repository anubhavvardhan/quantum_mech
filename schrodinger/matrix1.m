
%> The program <matrix1.m> solves for approximate eigenvalues and eigenfunctions 
%> using the matrix diagonalization methods built into MATLAB.
%> Hamiltonian = kinetic part + potential defined in [0,1] +
%> BC: u(0) = u(1) = 0.
%> The potential is chosen from a set defined in <pot1.m> using <choice.m>, 
%> and the strength of it is defined by a multiplier. The multiplier must be 
%> chosen large enough in order to get interesting effects. 
%> Recall the scale of the spectrum of a particle confined to [0,1] is 
%>  E = (n*pi)^2/2,  approx = 5, 20, 44, 79, 123, 178..
%> The basis set are harmonic functions in [0,1] adapted to the BCs, you
%> choose the number of them used, say from 50 up to 150 depending on your PC. 
%> The matrix elements of the kinetic energy are defined analytically, 
%> those of the potential are calculated using the FFT built into MATLAB.
%> Then the eigenvalues and eigenfunctions of the resulting finite matrix are
%> found using the MATLAB <eig> algorithm and displayed graphically.
%> POSSIBLE ERRORS: Too few  basis functions will result in easily observed   
%> errors in the order of the eigenfunctions: the well-known relation between 
%> the number of nodes and the order of the eigenvalue can then fail.
%> The so calculated eigenfunctions can then be used to find the evolution of
%> an initially Gaussian wave packet.
%> The resulting evolution will be interesting
%> if the potential and velocity are suitably chosen.
%> ================================
 %>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 971106

close;format long;  % global Q1 , global PP

disp('> Welcome to <matrix1>');

txt={' MATRIX APPROXIMATION FOR EIGENVALUES AND EIGENFUNCTIONS' 
' ' 
' This program solves for approximate eigenvalues and eigenfunctions' 
' for 1D motion in an interval [0,1] with a choice of potential.' 
' It uses the matrix diagonalization methods built into MATLAB.' 
' ' 
' The Hamiltonian = kinetic part + potential defined in [0,1]' 
' Boundary conditions u(0) = u(1) = 0.' 
' The potential is chosen from a set defined in the file <pot1.m>.' 
' The strength of it is defined by a multiplier.' 
' The multiplier must be chosen large enough in order to get' 
' interesting effects. ' 
' ' 
' See more information using the help in the file <matrix1>.'};

tt1=text0([.15 .4 .75 .45], txt);

cbutt(.75,.01);
delete(tt1);

Q2=2;

while  Q2==1|Q2==2

			if Q2==2

choice;  %% << Choice of potential in <pot1.m>.

P='pot1(X)';

G=200; % no of points in display of potential and wave functions
X=linspace(0,1,G); V=eval(P);
		
%%% calculating the FFT of the potential 
N=2^12; % no of points in FFT
X1=linspace(0,1-1/N,N); 
P1=strrep(P,'X','X1'); V1=eval(P1);		
X2=[X1(2:N),1];
P2=strrep(P,'X','X2');
V2=eval(P2);
VV=[V1,fliplr(V2)];
Y=fft(VV)/N/2; % cosine FFT 
Y1=real(Y);
%%%

	 	end % of Q2=2

Q3=2;

while Q3==2

%% Display of of potential:
subplot(2,2,3), l1=plot(X,V); set(l1,'LineWidth',3); 
title('This is your potential, normalized');

txt={' You must choose the number of harmonic functions used in the' 
' basis to construct the matrix approximation. ' 
' The calculations get slower when this number is larger than' 
' about 100. A low value gives bad accuracy.' 
' Now insert your choice!'};

tt1=text0([.15 .55 .75 .25], txt);

gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='50';ee=edit1([.75 .55 .15  .05],str);
N0=eval(ee); delete(tt1),

%%% calculate the matrix representations of potential and kinetic energy

		NN=1:N0; YY=Y1(NN);		
		pot = toeplitz(Y1(1:N0)) - hankel(Y1(3:N0+2));
		kin=pi^2 *NN.^2/2; kin=diag(kin);

txt={' You must choose the multiplier of your chosen potential. ' 
' Recall that the eigenvalues of the 1D box is approximately' 
'  5, 20, 44, 79, 123, 178.. ' 
' You had better choose a multiplier in the range 10^3 - 10^5' 
' in order to see interesting effects.' 
' Now insert your choice, including the SIGN!'};

tt1=text0([.15 .55 .75 .25], txt);

gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='1000'; ee=edit1([.75 .55 .15  .05],str);
mult=eval(ee); delete(tt1),

subplot(2,2,3);
l1=plot(X,mult*V); set(l1,'LineWidth',3); 
title('This is your potential with multiplier');

txt={' You must choose the number of eigenvalues/eigenfunctions' 
' you want to see displayed. ' 
' ' 
' This number should be much smaller than the number of modes' 
' chosen for the expansion.' 
' ' 
' Now insert your choice'};

tt1=text0([.15 .55 .75 .25], txt);

gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='10'; ee=edit1([.75 .55 .15  .05],str);
q5=eval(ee); delete(tt1),

%%% the matrix representation of total Hamiltonian

W=kin+mult*pot; % Hamiltonian as a N0 x N0 matrix

[K,D]=eig(W); 
%  [K,D] = eig(W) --> diagonal matrix D of eigenvalues and a
%  matrix K whose columns are the corresponding eigenvectors such
%  that W*K = K*D.
ww=real(D); % to remove spurious small imaginary parts
[ww,ID]=sort(diag(ww));
% [Y,I] = SORT(X) also returns an index matrix I. 
%  If X is a vector, then Y = X(I). 
 
ww1=ww(1:q5); % A subset of eigenvalues for display
disp(ww1);

% Bar diagrams over eigenvalues and difference vector
subplot(1,2,1), bar(ww1),
title(sprintf('Eigenvalues up to order #%g ',q5));
subplot(1,2,2), bar(diff(ww1)),
title('Differences between eigenvalues');
xlabel('Press the button to continue!');
cbutt(.75,.01);

% Energy levels superimposed on potential % figure(gcf);
subplot(1,1,1),
l1=plot(X,mult*V,'-.b');set(l1,'LineWidth',3);hold on; 
plot(X,ww1*ones(size(X)),'b'),hold off,
title('Lowest energy levels superimposed on potential');

YY=sqrt(2)*sin(pi*NN'*X); YY=K(:,ID)'*YY;

rbutt([0.55  0.01  0.2  0.05],'New parameters','uiresume; Q3=2;'),
gbutt([0.75  0.01  0.15  0.05],'CONTINUE','uiresume;Q3=3;'), 
uiwait; 
obutt;

end %of Q3

rbutt([.75 .01 .15 .05],'WAIT..',''), drawnow;

		for n=1:q5

		l1=plot(X,[YY(n,:);zeros(size(X))]);
		title(sprintf('Eigenfunction #%g ',n)),
		xlabel(sprintf('Eigenvalue = %g',ww(n))),
		set(l1,'LineWidth',2),
		drawnow,
		pause(1);
	
		end
		
cbutt(.75,.01);
xlabel('Press the button to continue!');

plot(X,YY(1:q5,:));
title(sprintf('All the %g eigenfunctions',q5));

rbutt([0.6  0.01  0.15  0.05],'Wave packet','uiresume; Q2=2;'); 
gbutt([0.75  0.01  0.15  0.05],'CONTINUE','uiresume; Q2=1;');
uiwait;

if Q2==2 % Show evolution of Gaussian wave packet


Q3=1;

while Q3==1

obutt,

L=1; % Length of interval

%%%%% INPUTS: center and velocity of wave packet
tt1=text0([.15 .5 .75 .1],'Enter the center of wave packet, in [0,1]');
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='.8'; ee=edit1([.75 .5  .15  .05],str);
delete(tt1), x0=eval(ee); 

txt=sprintf('You should choose a velocity smaller than %g ',N0/L);
tt1=text0([.15 .5 .75 .1],txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='100'; ee=edit1([.75 .5  .15  .05],str);
delete(tt1), k0=eval(ee); 

txt={'Enter the number of time steps in evolution'};
tt1=text0([.15 .5 .75 .1],txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='100'; ee=edit1([.75 .5  .15  .05],str);
delete(tt1), M=eval(ee);


sigma0 = L/30;   % Standard deviation of the wavefunction - change if you need
Norm = 1/(sqrt(sigma0*sqrt(pi)));  % Normalization
Y = Norm * exp(i*k0*X) .*exp(-(X-x0).^2/(2*sigma0^2)); %row vector
Y = Norm * exp(i*k0*X1) .*exp(-(X1-x0).^2/(2*sigma0^2));

Y(1)=0;Y=[Y,0,-fliplr(Y(2:N))];% Y(N)=0;
FFY=i*fft(Y)/N/2; % sine FFT  %FY(1)=0; %????????
FY=FFY(2:N0+1);FY=FY*K(:,ID);

% Parameters of time evolution
%tau=input('> Enter the time step per picture (start with 0.001) >>>');
tau=0.015/(abs(k0) + 20); 
%50 pictures for each crossing, discounting potential energy.

disp('> Calculating the evolved probability distributions, wait ...');
TT=tau*[0:M]; %row exp(i*diag(D)*TT).vector
FYT=((FY'*ones(size(TT))).*exp(i*ww*TT))'*YY;
Y7=abs(conj(FYT).*FYT);  % evolved probability density
q1=M;
norm=1/sum(Y7(1,:));
Y7=norm*Y7;

l1=plot(X,Y7(1,:));title('Initial distribution');xy=axis;
set(l1,'LineWidth',2);
xlabel('Press the button to continue!');
cbutt(.75,.01);
rbutt([.75 .01 .15 .05],'WAIT..',''), drawnow; 

for m=1:q1

l1=plot(X,Y7(m,:)); axis(xy);axis('off');	
title(sprintf(' Picture number %g ',m));
set(l1,'LineWidth',2);
drawnow; % 	pause(1)

end

xlabel('Press the button to continue!');
cbutt(.75,.01);
MX=Y7*X'; %average <X> as a function of time
MXX=sqrt(Y7*(X.^2)' - MX.^2); %variance - changed normalization.???

l1=plot(TT,[MX,MXX,1-sum(Y7')']);
set(l1,'LineWidth',2);
legend('Mean value','Deviation','1 - Norm');
title('Average and variance of position as functions of time ');

rbutt([.6  .01  .15  .05],'New packet','uiresume; Q3=1;');
gbutt([.75 .01  .15  .05],'CONTINUE','uiresume; Q3=2;'); uiwait;
 
end % of Q3

end % of Q2

rbutt([.2  .01  .2   .05],'New parameters','uiresume, obutt, Q2=1;');
rbutt([.4  .01  .2   .05],'New potential','close; Q2=2;');
bbutt([.6  .01  .15  .05],'BACK','close;Q2=3;'); 
bbutt([.75 .01  .15  .05],'QUIT','close;Q2=4;');
uiwait;

Q1=1;

	if Q2==3
	
	inbox; return;

	elseif Q2==2
	
	clear; Q2=2;

   end 

end % of Q2

 close;
	disp('> Type <matrix1> to do this again. Bye!'); 

%%% © Göran Lindblad 1997
