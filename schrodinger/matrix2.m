%> The program <matrix2.m> solves for approximate eigenvalues 
%> using the matrix diagonalization methods built into MATLAB.
%> Hamiltonian = kinetic part + potential defined in [0,1],
%> Boundary conditions: u(0) = u(1) = 0. The potential is chosen 
%> from a set defined in <pot1.m> using <choice.m>, and the 
%> strength of it is defined by a multiplier. The multiplier must  
%> be chosen large enough in order to get interesting effects. 
%> Recall the scale of the spectrum of a particle confined to 
%> [0,1] is  E = (n*pi)^2/2, approx = 5, 20, 44, 79, 123, 178, ..
%> The basis set are harmonic functions in [0,1] adapted to the BCs.
%> The matrix elements of the kinetic energy are diagonal and exact, 
%> those of the potential are calculated using the FFT. Then the 
%> eigenvalues the resulting finite matrix are ound using the  
%> MATLAB <eig> algorithm and a few are displayed graphically.
%> 
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 971106

close; format long; disp('> Welcome to <matrix2>');

txt={' MATRIX APPROXIMATION FOR EIGENVALUES' 

' ' 
' This program solves for approximate eigenvalues for a 1D ' 
' Schrodinger particle in an interval [0,1] with a choice of potential.' 
' It uses the matrix diagonalization routine of MATLAB.' 
' ' 
' The Hamiltonian = kinetic part + potential defined in [0,1]' 
' The boundary conditions are u(0) = u(1) = 0.' 
' The potential is chosen from a set defined in the file <pot1.m>.' 
' The strength of it is defined by a multiplier which can be changed' 
' using a slider.' 
' ' 
' For more information use the help in the file <matrix2>.'};

tt1=text0([.15 .4 .75 .45], txt);

cbutt(.75,.01); delete(tt1),

Q2=2;

while  Q2==1|Q2==2

	if Q2==2

choice;  %% << Choice of potential in <pot1.m>.

%%% FIXED PARAMETERS, CHANGE HERE IF YOU LIKE
G=200; % no of points in display of potential 
N0=150;% dimension of the matrix approx
q5=21; % number of eigenvalues displayed

%%% FFT of potential
P='pot1(X)';			
N=2^12; % no of points in FFT
X1=linspace(0,1-1/N,N); 
P1=strrep(P,'X','X1'); V1=eval(P1);		
X2=[X1(2:N),1];
P2=strrep(P,'X','X2');
V2=eval(P2);
VV=[V1,fliplr(V2)];
Y=fft(VV)/N/2; % this is the cosine FFT 
Y1=real(Y);

%% The matrix elements of the potential
NN=1:N0; %YY=Y1(NN);		
pot = toeplitz(Y1(1:N0)) - hankel(Y1(3:N0+2));

% Display of potential
X=linspace(0,1,G); V=eval(P); 
subplot(2,2,3), l1=plot(X,V); set(l1,'LineWidth',3); 
title('This is your potential, normalized');
% Choice of multiplier


txt={'You must choose the multiplier of the potential, in particular the SIGN!'};

tt1=text0([.15 .5 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='1000';	ee=edit1([.75 .5 .15  .05],str);
mult=eval(ee); delete(tt1),


%%% Now the kinetic part of the Hamiltonian
kin=pi^2 *NN.^2/2;	kin=diag(kin);

Q2=1;

	end % of Q2=2
		
vmin=min(mult*V);		
string0='The multiplier =  x ';
string=sprintf(strrep(string0,'x','%g'),mult); %replace string by value 
disp(string);

W=kin+mult*pot; % Hamiltonian as a N0 x N0 matrix

ww=eig(W); ww=real(ww);ww=sort(ww); % find and sort eigenvalues
 
ww1=ww(1:q5); % A subset of eigenvalues for display
% ww1=[0;ww1]; %%%% >????
disp('> The lowest eigenvalues are:');
disp(ww1);
disp('> The differences are:');
disp(diff(ww1));

% Displaying the eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1), bar(ww1); axis('tight');
title(sprintf('Eigenvalues up to order #%g ',q5));
legend(string);
subplot(1,2,2),
plot(X,[vmin;ww1]*ones(size(X)),'b'), axis('tight'),
hold on,
l1=plot(X,mult*V,'b');set(l1,'LineWidth',3), 
hold off,
title('Lowest energy levels superimposed on potential'),
xlabel('Press the green button to continue!');
cbutt(.75,.01);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1),
title('Choose multiplier using the slider!')
subplot(1,2,2),
bar(diff(ww1)),axis('tight');
title('Differences between eigenvalues');

rbutt([.55  .01  .2 .05],'New potential','close; Q2=2;');
gbutt([.75  .01  .15  .05],'CONTINUE','uiresume, Q2=3;');
g0 = slider0([.48  .12  .05  .8], [-2,2,0,], 'Q2=1;uiresume;');
uiwait; 

if Q2==1
mx=get(g0,'Value'); mult=mult*exp(mx), delete(g0), obutt, drawnow;
elseif Q2==3
delete(g0),
end


end %%%of Q2==1 or Q2==2

txt='> Type <matrix2> to do this again!';

rbutt([.3  .01  .15 .05],'REPEAT','close; matrix2; return;');
bbutt([.45   .01  .15  .05],'BACK','close; inbox;return;');
bbutt([.6   .01  .15  .05],'MAIN MENU','close; start; return;');
bbutt([.75  .01  .15  .05],'QUIT','close;disp(txt)');

