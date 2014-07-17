
%> The file <bound3d.m> finds the bound states for a central potential for each
%> value of the angular momentum L.
%> The potential is defined on [0,infinity] by an editable string.
%> The L-dependent 'centrifugal' part is added. 
%> The interval should be chosen such that the weakest bound states have
%> decayed to insignificant values at outside. 
%> The Fourier Grid Hamiltonian method is used.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 990414

clear; close; disp('> Welcome to <bound3d>');

q1=2; exx=[];

while q1 > 0
		
		if q1==2


txt= {' BOUND STATES FOR A CENTRAL POTENTIAL' 
' ' 
' The program <bound3d.m> finds the bound states for a central' 
' potential for each value of the angular momentum L. ' 
' The potential is defined on [0,infinity] by a string which can be'
' edited. The L-dependent centrifugal part is added.'
' In order to get correct values for weakly bound states'
' it is important to set the outer radius sufficiently large.'
' You can change this and other choices by editing the file.' 
'' 
' The Fourier Grid Hamiltonian method is used to construct a discrete'
' approximation on a finite interval for each value of L. ' 
' The program finds the bound state energies and stops when there are'
' no more bound states. '
''
' You can make the following choices:'
''
' (1) The spherical square well.'
' (2) The spherical harmonic oscillator (exactly solvable) '
' (3) Woods-Saxon potential '
' (4) Woods-Saxon and barrier '
' (5) A spherical square well and barrier'
' (6) A Yukawa-like potential '
''
' Choose by pressing the corresponding button!'
' ' 
' '};

N0 = 150; % the number of lattice points in the FGH method.
N1=N0+1;

t1= text0([.15 .2 .75 .75],txt);

gbutt([.3  .01 .1 .05],'#1','close;q0=1;');
gbutt([.4  .01 .1 .05],'#2','close;q0=2;');
gbutt([.5  .01 .1 .05],'#3','close;q0=3;');
gbutt([.6  .01 .1 .05],'#4','close;q0=4;');
gbutt([.7  .01 .1 .05],'#5','close;q0=5;');
gbutt([.8  .01 .1 .05],'#6','close;q0=6;');
uiwait, 



if q0==1
%% Spherical square well
x1=0; x2=12;
pot='1.5*(sign(X-4) -1) ';

elseif q0==2

%% Spherical harmonic oscillator
x1=0; x2=12;
pot=' .03*(X.^2 - 100)';

elseif q0==3

% Woods-Saxon potential 
x1=0; x2=12;
pot = '- 4./(1 + exp((X-3)/0.6))';  

elseif q0==4

% Woods-Saxon + Barrier
x1=0; x2=12;
pot= '4*(- 1./(1 + 1*exp((X-3)/.6)) + 1.5* exp(-.15*(X-5).^2))';


elseif q0==5
% spherical square well and square barrier 
x1=0; x2=12;
pot= '1.5*((sign(X-4) - 1)+ (sign(X-4) - sign(X-6)))'; 

elseif q0==6

% Yukawa-like potential 
x1=0; x2=12;

pot= '- 5*exp(- .2*X) ./(X + 1)';

end  

q1=1;

end 

if q1==1

X=linspace(x1,x2,N1);
xx=X;
% evaluate potential;
yy=eval(pot); 
miny=min(yy);


if miny >= 0

disp('> V(X) should have some negative part!');

bound3d; 

return

end 

w=[x1,x2,min(yy)-0.1,1.5*max(yy)+ 0.2*abs(min(yy))+0.1];

l1= plot(xx,yy);	axis(w);set(l1,'LineWidth',3);
title('This is the potential defined for L = 0');
xlabel('Press button to continue'); hold off;
cbutt(.75,.01);

rbutt([.75  .01  .15  .05],'WAIT...',''); drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% defining the potential string P

P='(V + 0.5*LL./((X+1.e-3).^2))';

P1=strrep(P,'V',pot); %% replacing V by string for potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=0;

zz=[];%% We start the iteration in L

	while L >= 0 
	
	% disp(sprintf('> L = %g',L));	
	
	P2=sprintf(strrep(P1,'LL','%g'),L*(L+1)); %% eff potential with L-value
	
	%% guessing the number of bound states
	% X=x;
	
	yy=eval(P2); %% value of potential
	
	w=wkb(xx,yy,.2,300);
	
	NB = length(w);
	
	% disp(NB)
	
	ee = fgh(P2,[x1,x2],N0,NB); 
	
	ind=find(ee <= 0);
	
	ee = ee(ind);
	
	z = ee';

	% disp(z); 
	
	ind = sum(~z);
	
	if ind == length(z) 
	
	z=[];
	
	else
		
	z=z(1:length(z)-ind); %% removes spurious zero solutions
		
	end
			
	if length(z)==0 %% no more roots stops the cycle
				
			L=-1;
				
	elseif length(z)==1
				
			if z <=  min(yy) + eps %% removes spurious solutions!
				
			L=-1;  %% no more roots stops the cycle
				
			else
					
			zz(L+1,2)=z;
									
			zz(L+1,1)=L;
			
			L=L+1;
					
			end

	else  %% z of lenght 2 or more

			if z(1) <= min(yy) + eps %% removes spurious solution!
		
			z=z(2:length(z));
					
			end
		
			zz(L+1,2:length(z)+1)=z;
				
			zz(L+1,1)=L;
				
			L=L+1;

	end %%  of if clause on length(z)

		end % end of L iteration
		
	[a,b]=size(zz);

			if a+b==0

			disp('> There are no bound states!');

			else

string2=str2mat('> RESULTS: ',...
'> The energies of the bound states are given below ',...  
sprintf('> for each value of the angular momentum L =0,..,%g',a-1),...
'');

null=~zz(:,2:b); null=~null;

minzz = floor(min(min(zz))-0.5);

[n1,n2]=size(null);

		if n2==1

		ind = null';

		else

		ind=sum(null'); %% indicates the number of bound states found  for each L

		end %% of if n2 

% delete(t1);

xy=[-0.2,length(ind), minzz ,1.1-0.02*minzz];
xx0=linspace(0,length(ind));
plot(xx0,0*xx0); drawnow,
title('Energy levels as functions of angular momentum L');
hold on;

		for M=1:length(ind)

		xx=0.8*linspace(0,1);
		xx1=ones(size(xx));		
		yyy=zz(M,[2:ind(M)+1])'*xx1;		
		s=plot(M+xx-1,yyy,'r');
		set(s,'LineWidth',2);
		axis(xy); 
text(M+0.1-1,0.5-0.01*minzz,sprintf('L=%g',M-1),'FontName','monaco','FontSize',12);
		%axis('off');

		end %% of M cycle
		
disp(string2);	disp(zz);

		end 
hold off;	

rbutt([.30  .01  .15  .05],'EDIT V(X)','uiresume, obutt, q1=1;');
rbutt([.45  .01  .15  .05],'REPEAT','uiresume, close, q1=2;');
bbutt([.60   .01  .15  .05],'BACK','uiresume, q1=3;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=0;'); 
uiwait; exx=[];

%rbutt([.45  .01  .15  .05],'REPEAT','uiresume; q1=1;'); 
%bbutt([.6   .01  .15  .05],'BACK','close; q1=2;');
%bbutt([.75  .01  .15  .05],'QUIT','close;q1=3;'); 
% uiwait;

	if q1==3

	boundst; return;
	
	elseif q1==1

% x1=-5; x2=5;
str2 = pot; 
gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume');
tt3=text0([.15 .2 .2 .05],'Edit the string >>>  ' );

str3=edit1([.35  .2  .55  .05],str2);

pot=str3;
delete(tt3),
	
	end	

	end
	
end

disp('> Type <bound3d> to do this example again!');

%%% © Goran Lindblad 1996
