
%> The file <special.m> provides a menu for the display of some of the most
%> frequently used special function, in particular some which are not in 
%> thee standard MATLAB toolbox = <specfun> directory/folder. 
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear; close;  disp('> Welcome to <special>!');

q=0;
 
axis('off')
axis([0 1 0 1])
t1=title('Special functions and expansions');
set(t1,'FontSize',18);
hold on

txt ={' '
	' These programs display the properties of some of the special      ' 
   	' functions used in solving the Schroedinger equation, in particular ' 
	' the Legendre and Bessel functions, and the harmonic oscillator     ' 
	' eigenfunctions. Expansions in these functions are also calculated. ' 
	' The standard MATLAB algorithms are used in most places, but        ' 
	' there are some exceptions.    '};

tt1=text0([.15 .45 .75 .4 ],txt); 

%%%%%%% INPUT

rbutt([.15  .36 .35  .06],'Legendre functions','close, q=1;')
rbutt([.15  .29  .35  .06],'Bessel functions','close, q=2;')
rbutt([.15  .22  .35  .06],'Spherical Bessels ','close, q=3;')
rbutt([.55  .36  .35  .06],'Expansion in Bessels','close,q=4;')
rbutt([.55  .29  .35  .06],'Harmonic oscillator','close,q=5;')
rbutt([.55  .22  .35  .06],'H-atom eigenfunctions','close,q=6;')
bbutt([.55  .15   .17   .06],'MAIN MENU','close,q=7;')
bbutt([.73  .15   .17   .06],'QUIT','close,q=8;')
uiwait;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if q==7
	start; return
	
	elseif q==8
	
	return
	
	elseif q==1 % Legendre polynomials

	legexp;
	
	elseif q==4 % Expansion in Bessels
	
	besjexp; return
	
	elseif q==5 % Harmonic oscillator eigenfunctions
	
	qho; return % 

	elseif q==6 % H-atom
	
	hatom;return	

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	elseif q==2 % Bessel functions
	

txt={'BESSEL FUNCTIONS'
''
'This file displays a finite number of Bessel functions'
'of integer order n and of the four different kinds '
'J_n , Y_n, I_n , K_n , as defined in HMF Chapter 9'
''
''
'Now choose an integer n >=  0 in the box' };


tt1=text0([.15 .5 .75 .25 ],txt); 
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='7';
ee=edit1([.75 .5 .15  .05],str);
N=eval(ee);

NN=[0:N-1]';
xmax=round(3+1.5*N);
x=linspace(0,xmax,200);
xy=[0,xmax,-1,1];
y=besselj(NN,x);
delete(tt1), plot(x,y) , axis(xy),
title(sprintf('The %g first Bessel functions J_n(x)',N));

y=bessely(NN,x);
cbutt(.75,.01);
plot(x,y); 	axis(xy); 
title(sprintf('The %g first Bessel functions Y_n(x)',N));
 
xy=[0,xmax/2,0,10];
y=besseli(NN,x);

cbutt(.75,.01);
plot(x,y);	axis(xy);
title(sprintf('The %g first Bessel functions I_n(x)',N));

y=besselk(NN,x);
cbutt(.75,.01);
plot(x,y); axis(xy);
title(sprintf('The %g first Bessel functions K_n(x)',N));

cbutt(.75,.01);

nn=linspace(0,10,50);
zz=linspace(0,25,50)';
ww1= besselj(nn,zz);

surf(nn,zz,ww1);
title('The family of Bessel functions J_{n}(x), for  n \in [0,10], x \in [0,25]');
ylabel('Coordinate x ');
xlabel('Parameter n ');
view(120,60);
cbutt(.75,.01);

ww= bessely(nn,zz);
ww1=crop(ww,1,-1);
figure(gcf);
surf(nn,zz,ww1);
title('The family of Bessel functions Y_{n}(x)  (values < -1 discarded) ');
ylabel('Coordinate x ');
xlabel('Parameter n ');
view(120,60);

cbutt(.75,.01);
special;	return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	elseif q==3 % Spherical Bessel functions
	
txt={'SPHERICAL BESSEL FUNCTIONS'
''
'The spherical Bessel functions can be calculated from the standard '
'ones in the standard MATLAB specfun directory routines. '
'The relation is given by HMT Chapter 10.1 :'
'       j(n,x) = sqrt(pi)*besselj(n+1/2,x)/sqrt(2*x), '
'and the same relation holds for the second kind functions.'
''
''
'Now insert an integer n „ 0 in the box' };

tt1=text0([.15 .5 .75 .3 ],txt); 
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='5';
ee=edit1([.75 .5 .15  .05],str);
N=eval(ee); delete(tt1),

xmax=10+N;
x=linspace(0, xmax)'+eps;
nn=[0:N-1];
y1=sqrt(pi)*besselj(nn+0.5,x);
y1=y1./(sqrt(2*x)*ones(1,N));
xy=[0,xmax,-0.5,1.1];
plot(x,y1);axis(xy);
title(sprintf('The %g first spherical Bessel functions j_n',N));
cbutt(.75,.01);

y2=sqrt(pi)*bessely(nn+0.5,x);
y2=y2./(sqrt(2*x)*ones(1,N));

xy=[0,xmax,-0.6,0.4];
plot(x,y2);axis(xy);
title(sprintf('The %g first spherical Bessel functions y_n',N));
cbutt(.75,.01);close;
special;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 	
	
	elseif q==0
	
	return;	
	
	end

	

