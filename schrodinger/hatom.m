
%> The file <hatom.m> gives a menu for programs displaysing some properties
%> of the  H atom eigenstates.
%>
%> © Goran Lindblad - gli@theophys.kth.se


close, clear, disp('> Welcome to <hatom>!');

q1=1;

axis('off');axis([0 1 0 1]), t1=title('States of the H atom');
set(t1,'FontSize',18);

ww = {''
	' This set of programs displays some properties of the bound states   ' 
	' of the  hydrogen atom, with graphics of both radial and angular' 
	' dependence.              ' 
	' ' 
	' The calculations are based on the MATLAB algorithms for the  ' 
	' Legendre functions and special algorithms for the radial parts.'};

text0([.15 .45 .75 .4], ww); 
rbutt([.15  .36  .35  .06],'Spherical harmonics','close,q=1;')
rbutt([.15  .29  .35  .06],'Angular dependence of orbitals','close,q=2;')
rbutt([.15  .22  .35  .06],'Radial wave functions','close,q=3;')
rbutt([.55  .36  .35  .06],'Atomic orbitals','close,q=4;')
rbutt([.55  .29  .35  .06],'Special functions','close,q=5;'),
bbutt([.55  .22  .17  .06],'MAIN MENU','close, q=6;'),
bbutt([.73  .22  .17  .06],'QUIT','close,q=7;'),
uiwait

if q==7
return;

elseif q==6
start; return;

elseif q==2 
angl;return;

elseif q==4
orbitals; return;

elseif q==5
special; return;

elseif q==1  %% Spherical harmonics
	

while q1==1

txt={' SPHERICAL HARMONICS DISPLAYED'
''
' The spherical harmonics are defined in terms of the associated '
' Legendre functions.'
''
' The algorithm calculating the Legendre functions is a standard '
' component of MATLAB.'
''
''
' Now insert a positive integer L in the box'};


tt1=text0([.15 .45 .75 .35], txt);

	
str='3';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .45 .15  .05],str);
L =eval(ee);delete(tt1);

%%%%%%% INPUT

x=linspace(0,1);
y1=ylm(L,x);
plot(x,y1);
xy=axis;
title(sprintf('The spherical harmonics of degree L = %g',L));
xlabel('The argument is cos(theta) in (0,1),  while phi = 0!');

str='The spherical harmonic Y(L,M)';

str1=sprintf(strrep(str,'L','%g'),L); 
rbutt([.75 .01 .15 .05], 'WAIT..', '');

pause(2) 

	for n=1:L+1

	l1=plot(x,y1(n,:)); 
	axis(xy);
	set(l1,'LineWidth',2); 
	str2=sprintf(strrep(str1,'M','%g'),n-1); 
	title(str2);
	xlabel('The argument is cos(theta) in (0,1),  while phi = 0!');

 	pause(1)
 
 	end 

plot(x,y1);
title(sprintf('The spherical harmonics of degree L = %g',L));
xlabel('The argument is cos(theta) in (0,1),  while phi = 0!');

cbutt(.75,.01);
disp('> We now display the same information in a polar diagram');
rbutt([.75 .01 .15 .05],'WAIT..',''), drawnow;
X=linspace(0,2*pi,200);

y1=ylm(L,cos(X));

str='The spherical harmonic Y(L,M), the amplitude squared is plotted ';

str1=sprintf(strrep(str,'L','%g'),L); 

for n=1:L+1

y2=y1(L+n,:).^2;

yy=y2.*cos(X); xx=y2.*sin(X);

l1= plot(xx,yy);axis('off');set(l1,'LineWidth',2);
axis('equal');xy=axis;

xmax=max(xy);
l1=line([0,0],1.1*[-xmax,xmax]);
set(l1,'LineWidth',3);
set(l1,'Color','r');
	str2=sprintf(strrep(str1,'M','%g'),n-1); 
	title(str2);

l1=text(0.1*xmax, 1.05*xmax, 'z-axis');

set(l1,'FontName','palatino'); set(l1,'FontSize',18);%set(l1,'Color','Black');
drawnow;

pause(1);

end % of iteration in n


rbutt([0.3 .01 .15 0.05],'REPEAT','close, q1=1;'); 
bbutt([0.45 .01,.15 0.05],'BACK','close, q1=2;');
bbutt([0.6 .01 .15 0.05],'MAIN MENU','close,q1=3;'); 
bbutt([0.75 .01 .15 0.05],'QUIT','close,q1=0;'); 
uiwait;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
		elseif q==3 %%%% Radial wave functions

q1=1;close;

while q1==1

txt={' RADIAL WAVE FUNCTIONS DISPLAYED' 
' ' 
' The radial wave functions are defined in terms of Laguerre ' 
' polynomials and are calculated using algorithms in <radial.m>' 
' and <laguerre.m>. .' 
' ' 
' ' 
' Now insert an integer n > 0 in the box'};


tt1=text0([.15 .45 .75 .4], txt);

	
str='3';
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
ee=edit1([.75 .45 .15  .05],str);
n =eval(ee);delete(tt1);


xmax=2+2*n*(n+1); x=linspace(0,xmax,300); y=radial(n,x');

plot(x,y);axis('tight'),
title(sprintf('Radial wave functions for states with n = %g', n));
cbutt(.75,.01);
y=(x'*ones(1,n)).*radial(n,x');
plot(x,y), axis('tight'),
title(sprintf('Radial wave functions with n = %g, multiplied by r', n)),
xlabel('Length scale in Bohr radi');

cbutt(.75,.01);
plot(x,y.^2), axis('tight'),
title(sprintf('Radial probability density for states with n = %g', n));
xlabel('Length scale in Bohr radi')
w=[];
		for k=1:n

		y=radial(k,x'); y= x'.*y(:,k); y=y.^2;	w=[w,y];

		end

cbutt(.75,.01);
plot(x,w);
title(sprintf('Radial probability density of the Rydberg (l=n-1) states up to n = %g', n));
xlabel('Length scale in Bohr radi');

rbutt([.3  .01  .15  .05],'REPEAT','uiresume; q1=1;'); 
bbutt([.45 .01  .15  .05],'BACK','uiresume; q1=2;');
bbutt([.6  .01  .15  .05],'MAIN MENU','close;q1=3;'); 
bbutt([.75 .01  .15  .05],'QUIT','close; q1=4;'); 
uiwait;

if q1 < 3
obutt,
end

end	

end	

	if  q1==3

	start;return 	

	elseif q1==2

	hatom; return
	
	end

disp('> Type <hatom> to repeat this display!');
