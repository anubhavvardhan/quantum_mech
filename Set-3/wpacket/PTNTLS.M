function [V, xo, Dx, Comment] = ptntls(x, potType)
%
%	function [V, xo, Dx, Comment] = PTNTLS(x, potType)
%
%	Generates variuos energy potential profiles in eV unit as
%	functions of x-axis positions and set wave-packet center and
%	x-indetermination to suitable values. The potential list can
%	be enriched by continuing the 'elseif' chain and introducing other
%	uimenu items in the routine WPACKET.CPP.  
%
%	INPUT:
%	x = coordinates of x-axis locations 
%	potType =  integer from 0 to 20 indicating a peculiar potential profile  
%		    (0 = flat; 1-5 positive boxes; 6-10 negative boxes;
%		    11 = positive linear; 12 = negative linear;
%		    13 = positive periodic; 14 = negative periodic; 
%		    15 = harmonic oscillator with canonical indetermination;
%		    16 = harmonic oscillator with sub-canonical indetermination;
%		    17 = harmonic oscillator with super-canonical indetermination;
%		    18 = symmetric quartic;
%		    19 = asymmetric quartic.)
%	 	
%	OUTPUT:
%	V = potential profile in eV
%	xo = location of wave-packet center at start
%	Dx = wave-packet indetermination 
%	Comment = string reporting the DX value (and its meaning) 
%
%				R.Nobili - Padova University, Dec. 99
	Vo=100;
	relMass =1;
	N=(length(x)-1)/2;
	pNum=2*N+1;

	if potType==0,
		V=zeros(size(x));
		xo=-10;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
	elseif potType==1,
		L=1;
		xo=-5;
		V=zeros(size(x));
		V((N-L):(N+L))=2*Vo*ones(1,2*L+1);
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];

	elseif potType==2,
		L=10;
		xo=-4;
		V=zeros(size(x));
		V((N-L):(N+L))=Vo*ones(1,2*L+1);
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];

	elseif potType==3,
		L=20;
		V=zeros(size(x));
		V((N-L):(N+L))=Vo*ones(1,2*L+1);	
		xo=-10;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];

	elseif potType==4,
		L=60;
		V=zeros(size(x));
		Vo=Vo;
		V((N-L):(N+L))=Vo*ones(1,2*L+1);	
		xo=-10;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
    
    elseif potType==5,
		L=80;
		V=zeros(size(x));
		Vo=Vo;
		V((N-L):(N+L))=Vo*ones(1,2*L+1);	
		xo=-10;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
        
	elseif potType==6,
		L=1;
		xo=-5;
		V=zeros(size(x));
		V((N-L):(N+L))=-2*Vo*ones(1,2*L+1);
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
        
	elseif potType==7,
		L=10;
		xo=-4;
		V=zeros(size(x));
		V((N-L):(N+L))=-Vo*ones(1,2*L+1);
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
        
	elseif potType==8,
		L=20;
		V=zeros(size(x));
		V((N-L):(N+L))=-Vo*ones(1,2*L+1);	
		xo=-10;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
        
	elseif potType==9,
		L=60;
		V=zeros(size(x));
		Vo=Vo;
		V((N-L):(N+L))=-Vo*ones(1,2*L+1);	
		xo=-10;
		Dx=0.3;
 	    Comment=['\Delta x = 0.3 ', 197];

	elseif potType==10,
		L=80;
		V=zeros(size(x));
		Vo=Vo;
		V((N-L):(N+L))=-Vo*ones(1,2*L+1);	
		xo=-10;
		Dx=0.3;
 	    Comment=['\Delta x = 0.3 ', 197]; 
	elseif potType==11,
		xo=-5;
		V=zeros(size(x));
		V((N+1):(2*N+1))=4*Vo*(1:N+1)/N;
        Dx=0.3;
 	    Comment=['\Delta x = 0.3 ', 197];

	elseif potType==12,
		xo=-5;
		V=zeros(size(x));
		V((N+1):(2*N+1))=-4*Vo*(1:N+1)/N;	
	elseif potType==13,
		L=3*ceil(N/4);
		V=zeros(size(x));
		Y=2*pi*2*(1:N+1-L)/(N-L);
		V(N+1:pNum-L)=Vo*sin(Y).^2;
		xo=-5;
		Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];
        
	elseif potType==14,
		L=3*ceil(N/4);
		V=zeros(size(x));
		Y=2*pi*2*(1:N+1-L)/(N-L);
		V((N+1):(2*N+1-L))=-Vo*sin(Y).^2;
		xo=-5;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];

	elseif potType==15,
		Kel=2; % elastic constant of harmonic oscillator in eV/Angstrom^2
		V = 0.5*Kel*x.*x; % harmonic oscillator potential
		h_bar = 0.6582122; % Plank's constant in eV*femtosec; 
		m_e=5.11e5; % electron mass in eV
		Mass=m_e*relMass;
		c = 2997.92458; % light velocity in Angstrom/femtosec
		% Dx=sqrt(c*h_bar/(2*sqrt(Mass*Kel))); %i.e. Dx=0.9879;
        Dx = 0.9864; % best behavior! Probably, routine discrete timing 
        % affects the canonical behavior with theoretical value of Dx. 
		% Thus, the wave function C*exp(-(x-xo)^2/4*Dx^2) has
		% the canonical indetermination for the QM harmonic oscillator    
		xo=0;
		Comment=['\Delta x = 0.9864 ', 197,' (canonical)'];
	elseif potType==16,
		Kel=2; % elastic constant in eV/Angstrom^2
		V = 0.5*Kel*x.*x; % harmonic oscillator potential
		Dx = 0.5; % sub-canonical indetermination
		xo=0;
		Comment=['\Delta x = 0.5 ', 197,  '(sub-canonical)'];
	elseif potType==17,
		Kel=2; % elastic constant in eV/Angstrom^2
		V = 0.5*Kel*x.*x; % harmonic oscillator potential
		Dx = 2; % super-canonical indetermination
		xo=0;
		Comment=['\Delta x = 2.0 ', 197,' (super-canonical)'];
	elseif potType==18,
		V = (x.*x.*x.*x)/200 - x.*x + 50; 
		Dx = 0.5; % 
		xo=-10;
		Comment=['\Delta x = 0.5 ', 197];
	elseif potType==19,
		V = (x.*x.*x.*x)/200 - x.*x -2.5*x + 50;
		Dx = 0.5; % 
		xo=-10;
		Comment=['\Delta x = 0.5 ', 197];
   elseif potType==20,
		V = (x.*x.*x.*x)/100 - 2*x.*x - 10*x + 10;
		Dx = 0.5; % 
		xo=-7.5;
		Comment=['\Delta x = 0.5 ', 197];
      
    elseif potType==21,
     	y=randn(size(x));
     	Nx=length(x);
        xc=x(floor(0.8*Nx));  
       Dxc = x(floor(Nx/20))-x(1);
       xo=-10;
       Dx= 0.5;
       G=exp(-(x-xc).*(x-xc)/(2*Dxc*Dxc));
       V=Vo*G.*y;
       Comment=['\Delta x = 0.5 ', 197];
    else
		V=zeros(size(x));
		xo=0;
        Dx=0.9;
 	    Comment=['\Delta x = 0.9 ', 197];

    end
		
