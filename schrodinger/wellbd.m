%> The file <wellbd.m> calculates bound states in spherical square 
%> potential well of finite depth.  The radius is chosen in CW, value A, thus:
%> V(r) = V < 0 for r < A, = 0 for r > A.  V is chosen in CW.
%> The bound state energy levels are found from the known analytical formulas 
%> involving Bessel functions by an automatic numerical search for the solutions 
%> of a transcendental equation. This search may miss weakly bound states.
%> The accuracy is largely determined by the parameter NN, the number of energy 
%> lattice points. The analytic formulas can be found e.g. in Messiah Sec. IX.10.
%> Uses standard Bessel function program of MATLAB, and <findzero.m>.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

clear; close; disp('> Welcome to <wellbd>!');

q1=1;

txt={' SPHERICAL SQUARE WELL: BOUND STATES' 
' ' 
' This file <wellbd.m> calculates bound states in spherical square'  
' well of finite depth.  The the depth and radius is chosen below. ' 
' ' 
' The bound state energy levels are found from the known analytical '  
' formulas involving Bessel functions. There is an automatic numerical' 
' search for the solutions of the boundary conditions.'}; 

tt0=text0([.15 .5 .75 .3], txt);
tt1=text0([.15 .45 .75 .05], 'Choose V negative at right:');
str='-3';
gbutt([.75  0.01  0.15  0.05],'CONTINUE','uiresume'); 
ee=edit1([.75 .45 .15  .05],str);
V=eval(ee);  delete(tt1);

tt2=text0([.15 .45 .75 .05], 'Choose radius A positive at right:');
str='4';
gbutt([.75  0.01  0.15  0.05],'CONTINUE','uiresume'); 
ee=edit1([.75 .45 .15  .05],str);
A=eval(ee);		
delete(tt0,tt2), 

while q1==1

rbutt([.75 .01 .15 .05], 'WAIT',''), drawnow,

% delete(tt1,tt2),

	if V > 0 

	disp('> The potential must be negative, changing sign!');	V=-V;

	end
		
	if A < 0 

	disp('> The radius must be positive, changing sign!');	A=-A;

	end
		
disp(sprintf('> V = %g',V));		
disp(sprintf('> A = %g',A));

% A rough estimate of maximal relevant value of L
Lmax=ceil(A*sqrt(-2*V))+1; 

disp(sprintf('> The maximal angular momentum Lmax = %g',Lmax));

NN=500; % no of points in energy lattice, CHANGE HERE IF YOU LIKE.

% We choose to make the lattice equispaced in the wave number of the
% (exponentially decaying) solutions outside the well:

k1=linspace(1,0,NN)*sqrt(-2*V+eps); %
unk=ones(size(k1));
k0=real(sqrt(-2*V - k1.^2 +eps));%% wave number inside well

LL0=[0:Lmax]+1; 

 LL=[0:Lmax+1]';

unl=ones(size(LL0))';


%% The regular radial solutions inside well are the spherical bessels
%% j_l(k0*r), the solutions regular outside are k_l(k1*r) 
%% The continuity of the logarthmic derivative at r = A reads 
%% as follows, with D as the derivative with respect to the 
%% variable u=k*r: k0*Dj_l / j_l = k1*Dk_l /k_l
%% The derivative follows from the recursion relations
%% D j_l = - j_(l+1) + l* j_l / u  
%% and we obtain the relation to solve for eigenvalues
%% k0 * j_(l+1) * k_l  - k1 * k_(l+1) * j_l


bess0=besselj(LL+1/2, A*k0)'; 

bess1=besselk(LL+1/2, A*k1)'; 

wro= (unl*k0).*bess0(LL0+1,:).*bess1(LL0,:)...
 - (unl*k1).*bess1(LL0+1,:).*bess0(LL0,:) ; 


% bess0=besselj(LL+1/2, A*k0)'; 
% nn1=(ones(size(LL))*k1*A).^((LL+1)*unk);
% bess1=nn1.*besselk(LL+1/2, A*k1)'; %% the regular solutions  outside 
% dbess0 = (LL1*unk).*bess0(1:Lmax+1,:) - (unl*k0*A).*bess0(2:Lmax+2,:);
% dbess1 =(LL1*unk).*bess1(1:Lmax+1,:) - bess1(2:Lmax+2,:);
% wro=bess0(1:Lmax+1,:).*dbess1 - bess1(1:Lmax+1,:).*dbess0;

zmax=ceil(sqrt(-2*V)/pi) ;

z=[]; n=1;

	while n > 0

	z=findzero(k1,wro(n,:)); 	z=-z.^2/2; 		ind = sum(~z); 
		
		if ind == length(z) 
		
		z=[];

		else
		
		z=z(1:length(z)-ind); %% removes spurious zero solutions
		
		end
		
		if length(z)==0 %% no more roots stops the cycle
				
		n=0;
				
		elseif length(z)==1
				
			if z <= V + eps %% removes spurious solutions!
				
			n=0;  %% no more roots stops the cycle
				
			else
					
			zz(n,2)=z;  zz(n,1)=n-1; 	n=n+1;
					
			end
						
				
		else  %% z of lenght 2 or more
				

			if z(1) <= V + eps %% removes spurious solution!
		
			z=z(2:length(z));
					
			end
		
		zz(n,2:length(z)+1)=z;		zz(n,1)=n-1;		n=n+1;
			
		end %%  of if cycle

	end  %% of n cycle


	w=size(zz); w=w(1)-1;

	if w < 0

	disp('> There are no bound states!');

	else		

	string2=str2mat('> RESULTS',...
	'> The energies of the bound states are given below:  ',...
	sprintf('> for each value of the angular momentum L =0,..,%g',w),...
	' ');

	disp(string2);	disp(zz);

	[a,b]=size(zz);

	null=~zz(:,2:b); null=~null;

	minzz = floor(min(min(zz)));

	[n1,n2]=size(null);

		if n2==1

		ind = null';

		else

		ind=sum(null'); %% indicates the number of bound states found  for each L

		end %% of if n2 


	xy=[-0.2,length(ind), minzz - 1,1-0.02*minzz];
	xx0=linspace(0,length(ind));
	plot(xx0,0*xx0);
	hold on;

		for M=1:length(ind)

		xx=.8*linspace(0,1);		xx1=ones(size(xx));
		
		yy=zz(M,[2:ind(M)+1])'*xx1;
		
		s=plot(M+xx-1,yy,'r');
		set(s,'LineWidth',2), title('Bound states as functions of the angular momentum L'),
		text(M+0.1-1,  0.5,sprintf('L=%g',M-1),'FontName','monaco','FontSize',12);
		axis(xy); 		%axis('off');

		end %% of M cycle

	end %% of if w statement
		
hold off;

rbutt([.3   .01  .15  .05],'REPEAT','uiresume, q1=1;'); 
bbutt([.45  .01  .15  .05],'BACK','uiresume, q1=2;');
bbutt([.6   .01  .15  .05],'MAIN MENU','uiresume,q1=3;'); 
bbutt([.75  .01  .15  .05],'QUIT','close,q1=0;'); 
uiwait;
	
	if q1==1
		
	obutt,	
	tt1=text0([.15 .15 .75 .05], 'Choose V negative at right:');
	str='-5';
	gbutt([.75  .01  .15  .05],'CONTINUE','uiresume'); 
	ee=edit1([.75 .15 .15  .05],str);
	V=eval(ee); delete(tt1);

	tt2=text0([.15 .15 .75 .05], 'Choose radius A positive at right:');
	str='2';
	gbutt([.75  .01  .15  .05],'CONTINUE','uiresume'); 
	ee=edit1([.75 .15 .15  .05],str);
	A=eval(ee);delete(tt2);

	elseif q1==2

	boundst; return;
		
	elseif q1==3
	
	start; return;
		
	end
		
end

close,

disp('> Type <wellbd> to do this problem all over again!');


