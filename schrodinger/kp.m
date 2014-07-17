
%> The file <kp.m> solves for eigenvalues of Kronig-Penney type periodic 1D potential.
%> It approximates the Hamiltonian by a finite matrix eigenvalue problem.
%> Periodic (Bloch) boundary conditions.  Potential chosen from list in <pot1.m> 
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

clear; close; disp('> Welcome to <kp>');

string={' KRONIG-PENNEY TYPE MODEL' 
' ' 
' This program solves for the lowest eigenvalues of a Kronig-Penney' 
' type model in 1D using a matrix approximation. ' 
' ' 
' There is a choice of potential in the unit cell of length 1.' 
' ' 
' The eigenvalues are displayed as wave numbers as function of the ' 
' Bloch wave number. Consequently, for the free particle you will' 
' get intersecting straight lines, not the parabolas you will ' 
' usually see in the textbooks.'};

tt1=text0([.15 .4 .75 .45], string);

cbutt(.75,.01); delete(tt1), 

Q=2;

		while Q > 0

		if Q==2

choice; % For choice of potential go to <choice.m>.

X0=linspace(0,1);
V0=pot1(X0);
subplot(2,2,3),
l1=plot(X0,V0);%only for display
set(l1,'LineWidth',3);
title('This is your potential, without multiplier!')
xlabel('Press the button to continue!');

%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%
txt={'Choose the size of Bloch momentum grid (10-50)'};
tt1=text0([.15 .5 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='15';
ee=edit1([.75 .5 .15  .05],str);
G2=eval(ee); delete(tt1),

%
txt={'Choose the number of eigenvalues you want to see '};
tt1=text0([.15 .5 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='7';
ee=edit1([.75 .5 .15  .05],str);
gg=eval(ee); delete(tt1),

bb=4*gg;%number of basis functions
G=4*bb;%sufficiently(?) large number of points to resolve basis functions
g2=ceil(sqrt(bb+1));% just keep a few  around the diagonal, there is a fast decay
% of matrix elements
GG=[0:g2];
GR=[g2+1:2*bb];
xx=linspace(0,1,G2+1); %Bloch grid
XX=linspace(0,1,G); %space grid
nn=-bb:bb;%grid of basis functions
yy=[exp(2*i*pi*GG'*XX);zeros(size(GR'))*XX];

		end 
%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%
txt={'Choose the multiplier of the potential, including SIGN!'};
tt1=text0([.15 .5 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='50';
ee=edit1([.75 .5 .15  .05],str);
multiplier=eval(ee);
subplot(2,2,3),
l1=plot(X0,multiplier*V0);%only for display
set(l1,'LineWidth',3);
title('This is your potential, with multiplier!'), delete(tt1),
rbutt([.75  .01  .15  .05],'WAIT..',''), drawnow;

V=multiplier*pot1(XX); %potential on space grid
Vmin=min(V); % minimum of potential 
vv=yy*V'/(G-1); % use FFT instead!!!!!!
VV=sparse(toeplitz(vv));
% VV=sparse(VV);

ww=[]; %start on matrix of eigenvalues 
w2=[]; %start on a vector of eigenvalues

%tic
for kappa=0:G2 %iterative solution over momentum grid

%disp(kappa);

k=pi*kappa/G2; % Bloch momentum

kk=2*pi*nn + k; 

kk=sparse(diag(kk.^2)./2);% kinetic energy matrix in basis grid

% kk=sparse(kk);

dd=eig(full(kk+VV)); %eigenvalues %%%  WE WOULD LIKE A REAL MATRIX

dd=real(dd);% make sure they are real

dd = sort(dd); %sorted eigenvalues

%nn=size(dd);

dd = dd(1:gg);

dd = real(sqrt(2*dd - Vmin)); % Wave numbers from bottom of band. 

ww=[ww,dd]; w2=[w2;dd]; w2=sort(w2); 

end

subplot(1,1,1),
plot(xx,ww');
title('Wave number eigenvalues (from bottom of band) as functions of Bloch wave number');
xlabel('Press the button to continue!');
cbutt(.75,.01);
		
stairs(w2);
title('Wave number eigenvalues (from bottom of band), ordered by value');
ylabel('Wave number');
xlabel('Order of eigenvalue. Large vertical steps correspond to forbidden bands');

rbutt([.2  .01  .2  .05],'New multiplier','uiresume, Q=1;'); 
rbutt([.4  .01  .2  .05],'New potential','close; Q=2;');
bbutt([.6  .01  .15 .05],'BACK','close;Q=3;'); 
bbutt([.75 .01  .15 .05],'QUIT','close;Q=0;'); 
uiwait; 
	
		if Q==3
		
		per1d; return;
		
		elseif Q==1
		
		obutt,
		
		end;		

end

disp('> Type <kp> to do this again. Bye!');



