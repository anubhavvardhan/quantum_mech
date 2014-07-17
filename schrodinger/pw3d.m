
%> The file <pw3d.m>  calculates partial, total, and differential cross
%> sections for a central potential specified by a string, and for
%> a range of wave numbers for the incident particle. 
%> You must choose a coordinate where the potential has reached it´s
%> asymptotic (zero), value, and where boundary conditions can be applied.
%> Also choose a range of wave numbers. 
%> It calculates the partial cross section for a given L, starting with
%> L = 0, you can step forward in L and stop when the cross sections become
%> negligible in the chosen interval.
%> The program can also calculate the first Born approximation.
%> The radial Schrodinger equation, with potential defined on [0,infinity]
%> by defined by a string, and the L-dependent centrifugal part added. 
%> BC: consistent with Bessels at x=0 and x=A. The parameter A should be 
%> chosen such that the potential is essentially zero outside x=A.
%> Integration method: Numerov algorithm.
%>
%> © Goran Lindblad - <gli@theophys.kth.se>

%> GL 961123

clear; close;q2=2; disp('> Welcome to <pw3d>!');

q1=1;
			
			while q1 > 0
			
			
			if q1==1

txt={' CENTRAL POTENTIAL: CROSS SECTIONS' 
' ' 
' The program <pw3d.m> calculates partial, total, and differential cross '
' sections for a central potential specified by a string and for a range'
' of wave  numbers. It calculates the partial wave shift for each value '
' of L, starting with L = 0, by numerical integration of the radial '
' equation using the Numerov method. It then finds the partial cross '
' section. You can step forward in L until the cross section becomes  ' 
' negligible in the chosen interval. The program can also calculate the'
' first Born approximation.' 
''
' You must choose the potential and a coordinate where the potential '
' has reached the asymptotic zero value; there the continuity condition'
' can be applied. You can also choose a range of wave numbers. ' 
' For more information, use the help for <pw3d>.'
''
' You can make the following choices:'
''
' (#1) The spherical square well.'
' (#2) Woods-Saxon potential '
' (#3) Woods-Saxon and barrier '
' (#4) A spherical square well and barrier'
' (#5) A Yukawa-like potential '
''
' Choose by pressing the corresponding button!'
' ' 
' '};


tt1=text0([.15 .2 .75 .75], txt);

%% INPUT = press one of the buttons

gbutt([.3  .01 .1 .05],'#1','uiresume;q0=1;');
gbutt([.4  .01 .1 .05],'#2','uiresume;q0=2;');
gbutt([.5  .01 .1 .05],'#3','uiresume;q0=3;');
gbutt([.6  .01 .1 .05],'#4','uiresume;q0=4;');
gbutt([.7  .01 .1 .05],'#5','uiresume;q0=5;');

uiwait, delete(tt1), obutt,

%% Defining each potential as a string <pot>

if q0==1
%% Spherical square well

pot='(sign(X-4) -1)';

elseif q0==2


% Woods-Saxon potential 

pot = '(- 2./(1 + exp((X-3)/.6)))';  

elseif q0==3


% Woods-Saxon + Barrier

pot= '(- 2./(1 + exp((X-3)/.6)) + .5* exp(-.25*(X-6).^2))';

elseif q0==4

% spherical square well and square barrier 

pot= '((sign(X-4) - 1)+ (sign(X-4) - sign(X-5)))'; 

elseif q0==5


% Yukawa-like potential 

pot= '2*(- exp(- .2*X) ./(X + .1))';

end  

str0=10; str1=1; str2='[0,2]';

q1=2;

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if q1==2

%% Display the potential defined by the the string. 

rr=150; rmax=10; 

X=linspace(0,rmax,rr);	Y=eval(pot);

w=[X(1),X(rr),min(Y)-0.1,1.5*max(Y)+ abs(min(Y))+0.1];
plot(X,Y);axis(w);
title('This is the potential for V = 1 and L=0');

%% INPUT = chooose value of x outside of which the potential is 
%% effectively = 0, using edit box

txt={'Choose the coordinate where the phase shift is determined >'};
% xlabel(txt);

tt1=text0([.15 .8 .75 .1], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');

ee=edit1([.75 .8 .15  .05],str0);
A=eval(ee); delete(tt1),
% disp(A);
q1=3;
	end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
	if q1==3

%% INPUT = scale the potential, using edit bxox

txt={'Choose a multiplier V for the potential (with sign!) > '};
tt1=text0([.15 .8 .75 .05], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str1=edit1([.75 .8 .15  .05],str1);
V=eval(str1); delete(tt1),
% disp(V);
			
	end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT = interval of wave numbers, you can choose this to study
%%  an interesting interval - edit box

txt={'Choose a wave number interval [k0,k1], (0 < k0 < k1 ) '};
tt1=text0([.15 .8 .75 .05], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str2=edit1([.75 .8 .15  .05],str2);
KK=eval(str2); delete(tt1),
% disp(KK);
rbutt([.75 .01  0.15  0.05],'WAIT..',''), drawnow;
% return 

yy=sum(Y)/(rr+0.5);
average=V*yy*rmax/A; %if this is negative we should decrease the lattice spacing. 
av=0.5*(average - abs(average)); kmax=KK(2); 
% Lmax=ceil(0.5*kmax*A+5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for the numerical integration

N=ceil(20*A*sqrt(kmax^2 - av)+100); % You can change here!
M=300; % The number of points in wave number range, change here!
			
cross=[]; ampl=[];

waveno= linspace(KK(1)+0.0001,KK(2),M); %% Lattice in wave number
K= 0.5*waveno.^2; %% Kinetic energy 
K1= ones(size(K)); 
K0= zeros(size(K));
delta=0.001*A; %% Starting point near singularity. You can change!
%%%%%%%%%% CHANGE HERE 
%%E0=V*pot2(delta); %% The (negative) value of pot at starting point;
E0=V*Y(1); %% The (negative) value of pot at starting point;

wave0=sqrt(2*(K-E0)); %% Local waveno at starting point
dd=wave0*delta;
h=A/(N-1); %% integration step
X=delta +[0:h:A]'; % range [delta, A + delta];
uniX=ones(size(X));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=A+delta; %% Right end point!
%%%%%%%%%%%% CHANGE HERE
%%%%%%%%% introducing the angular momentum part and scaling by V
P='V*P0+0.5*LL./((X+eps).^2)'; 
P1=strrep(P,'P0',pot); 
P2=sprintf(strrep(P1,'V','%g'),V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q3=1; L=0;

% 	for L=0:Lmax %% Calculation of cross section for given L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		while q3 > 0

P1=sprintf(strrep(P2,'LL','%g'),L*(L+1));
NN=waveno.^L;

	if L==0
	Y0=sin(dd); DY0=wave0.*cos(dd);
	NN=1;

	elseif L==1
	Y0=dd.*(sin(dd)-dd.*cos(dd));
	DY0=wave0.*(dd.*cos(dd) + (dd.^2-1).*sin(dd));
	
	else		
	Y0=delta*K1.*(1-dd.^2/2/(2*L+3)); 
	DY0=K1.*(L+1 - (L+3)*dd.^2/2/(2*L+3));
	NN=1./(1+1./NN);
	end

%% Integration
	ww=numerov1(delta,A,N,P1,1,K,Y0,DY0);
	
%% Free solutions in x =  A
j1=sqrt(A*waveno).*besselj(L+1/2,A*waveno);
y1=NN.*sqrt(A*waveno).*bessely(L+1/2,A*waveno);
%% Following functions are D (x j1(x)) etc for x = A*waveno
j2=(L+1)*j1/A - waveno.*sqrt(A*waveno).*besselj(L+3/2,A*waveno);
y2=(L+1)*y1/A - NN.*waveno.*sqrt(A*waveno).*bessely(L+3/2,A*waveno);

%% Find sin of phase shift
aa=y2.*ww(1,:) - y1.*ww(2,:);
bb=NN.*(j2.*ww(1,:) - j1.*ww(2,:));
normw=sqrt(aa.^2 + bb.^2 );
sinw=bb./normw; %% sin  of PW phase shift

phase=atan(bb./aa);phase=glue(phase,pi);

subplot(2,1,1),plot(waveno,phase);
axis tight, title(sprintf('Phase shift for L = %g.',L));
%% Cross section and amplitude
crs=4*pi*(2*L+1)*(sinw./waveno).^2; %% PW cross section
cross=[cross;crs]; %% Collect the PW cross sections
amp=(2*L +1).*sinw.*(aa + i*bb)./waveno./normw;
ampl=[ampl;amp];


%% Display the partial cross section % 
subplot(2,1,2), plot(waveno,crs),axis tight,
title(sprintf('Partial cross section for L = %g.',L)),
xlabel('Wave number');

gbutt([.60  .01  .15  .05],'L => L+1','uiresume;q3=1;'),
bbutt([.75  .01  .15  .05],'ENOUGH','uiresume;obutt;q3=0;'),
uiwait; 

if q3==1
L=L+1;
end

	end  %% of iteration in L governed by q3
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
% Plot the total and the partial cross sections
subplot(1,1,1);

if L==0
totcross=cross;
else
totcross=sum(cross);
end

l0=plot(waveno,[totcross;zeros(size(waveno))],'k');axis tight; 
set(l0,'Linewidth',2);
hold on
plot(waveno,cross);
axis tight;
title(sprintf('Total and partial cross sections up to Lmax = %g',L));
xlabel('Wave number'); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculating differential cross  sections
xlabel('Press the button to continue!');
cbutt(.75,.01);

theta=linspace(0,pi,50);
LL=[0:L];
pol=legpol(L,cos(theta));
ampl=pol'*ampl;
difcros=log(real(conj(ampl).*ampl)+1);
mesh(theta*180/pi,fliplr(waveno),rot90(difcros)); axis('tight'),
title('Diff. cross section as a function of theta and wave number, logarithmic scale!');
ylabel('Wave number'); xlabel('Angle theta');
view([1.3,1,5]);


rbutt([.3   .01  .15  .05],'Born appr','uiresume, obutt, q2=1;'); 
rbutt([.45  .01  .15  .05],'Set [k0,k1]','uiresume, obutt, q2=0;q1=4;'); 
rbutt([.6   .01  .15  .05],'Set V','uiresume, obutt, q2=0;q1=3;');
gbutt([.75  .01  .15  .05],'CONTINUE','uiresume, q2=2;q1=0;'); 
uiwait;

	if q2==1
	
	% Calculating the Born approximation
	disp('> Calculating the Born approximation');
	disp('> This may take some time!!');
	rbutt([.75 .01 .15 .05], 'WAIT...','');
	
	D=100; % This is a parameter you can change - smaller is faster!
	x=linspace(0,A,D);
	unix=ones(size(x));
	y=V*pot2(x); %%%%%%%%%%%%change here
	Lmax=L;
	L=0;

	tcross=zeros(size(waveno));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		for L=0:Lmax

		bess=besselj(L+1/2,x'*waveno).^2;
		phase=-A*pi*(y.*x)*bess/D;
		bcross=4*pi*(2*L+1)*(sin(phase)./waveno).^2;
      	tcross=tcross+bcross;
% plot a comparision of the calculated cross section and the Born approx.
		figure(gcf);
		plot(waveno,[cross(L+1,:);bcross])
		title(sprintf('Partial cross section and Born approximation for L=%g', L));
		xlabel('Wave number');	drawnow;
		end 		
	
		
cbutt(.75,.01);
plot(waveno,[totcross;tcross]), axis tight;
title(sprintf('Total cross section and Born approximation up to L = %g', Lmax));
xlabel('Wave number'),drawnow;

	q2=2; end %% q2 == 1

		if q2==2
		
% rbutt([.15  .01  .15  .05],'Set [k0,k1]','uiresume, obutt, q1=3;'); 
rbutt([.3   .01  .15  .05],'Edit V(X)','uiresume, obutt, q1=2; ');
rbutt([.45  .01  .15  .05],'REPEAT','close,  q1=1; ');
bbutt([.6   .01  .15  .05],'BACK','close, q2=3; q1=0;'); 
bbutt([.75  .01  .15  .05],'QUIT','close, q1=0;'); 
uiwait; % close;

		end % on q2 == 1 condition
		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		if q1==2
		
pp1 = pot; 

gbutt([.75 .01 .15 .05],'CONTINUE', 'uiresume');
tt3=text0([.15 .2 .2 .05],'Edit the string >>>  ' );

pot=edit1([.35  .2  .55  .05],pp1);

delete(tt3),	
		
		end %% of q1==1

		end  %%% of q1 condition

			if q2 == 3

			scatt3d;	return;			
						
			end 

disp('> Type <pw3d> to do this example again!');


