
%> The file <legexp.m> gives a demo of the Legendre polynomials and 
%> associated Legendre functions. 
%> It also calculates expansions of functions in such a ON basis.
%> The definition of the Legendre polynomials can be found in most 
%> textbooks on quantum mechanics or angular momentum.  
%> Here we use the files <legpol.m> and <legfun.m> for rapid calculation
%> of matrices of function values. 
%> MATLAB also has a standard algorithm for the Legendre functions, 
%> which includes the polynomials as a special case. 
%> Notation P(M,L), L >= M >= 0 integers. 
%> M = 0 gives the Legendre polynomials.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

clear; close; disp('> Welcome to <legexp>');


q=3; %%%%%%%%%%%%%%%

		while q > 0

		if q==3;%%%%%%%%%%%%%%%%

txt={'LEGENDRE POLYNOMIALS AND FUNCTIONS'
''
'This file gives a demo of the Legendre polynomials and associated'
'functions, and calculates expansions in series of such functions.'
'The definition of the Legendre polynomials can be found in most '
'textbooks on quantum mechanics or on angular momentum. '
'Here we use the files <legpol.m> and <legfun.m> for rapid '
'calculation of matrices of function values. '
''
'Notation P(M,L), L >=  M >=„ 0 integers. '
'M = 0 gives the Legendre polynomials.'};

tt1=text0([.15 .5 .75 .35], txt);

str='500'; %%%%%%%%%%%%%%%%%%%
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume');
tt2=text0([.15 .45 .75 .05],'Set the number of lattice points for integration');
ee=edit1([.75 .45 .15  .05],str);
nn=eval(ee); delete(tt2),
		
str='0'; %%%%%%%%%%%%%%%%%%%%
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume');
tt2=text0([.15 .45 .75 .05],'Set the index M of P(M,L)');
ee=edit1([.75 .45 .15  .05],str);
M=eval(ee); delete(tt1,tt2),
		
N=7; %%%%%%%%%%%%%%%%%% change here 
NN=1:N;
figure(gcf);
X=linspace(-1,1);

if M==0
w=legpol(N-1,X); 
string='The %g first Legendre polynomials';
string1='The Legendre polynomial of degree %g';

else

w=legfun(M+N-1,M,X);
string='The %g first Legendre functions';
string1='The Legendre function of degree %g';

end

plot(X,w); title(sprintf(string,N)); xy=axis;
cbutt(.75 ,.01);

q=2;

end % of q=3;

if q==2  %%%%%%%%%%%%%%%%%%%%%%%%%%

str=' abs(X+.5).*abs(.6-X + abs(.6 - X))';%%%%%%%%%%%%%%%%%%%%%
gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume');
tt2=text0([.15 .45 .75 .05],'Give a real function of X on [-1,1]');
ee=edit1([.15 .4 .75  .05],str);
string=ee;

delete(tt2),

X=linspace(-1,1,nn);X1=ones(size(X));

figure(gcf);y=eval(string);
ymax=max(y); ymin=min(y); delta=ymax-ymin;
xy=[-1, 1, ymin - 0.1*delta, ymax + 0.1*delta]; 
l1=plot(X,y,'r');
axis(xy);
set(l1,'LineWidth',2);
title('This is the function to be expanded');

N0=50;%%%%%%%%%%%%%% max no of terms in expansion, change here %%%%

NN=[1:N0]; yy=ones(size(NN'))*y;

if M==0 %%%%%%%%%% Polynomial case

w=legpol(N0-1,X); cc=yy.*w;
cc=trapz(X',cc'); % Integration by trapezoidal rule
cc=cc.*(2*NN-1)/2;% Normalization

else  %%%%%%%%%%% general case

w=legfun(M+N0-1,M,X);cc=yy.*w;
cc=trapz(X',cc'); % Integration by trapezoidal rule

end

cc1=cc'*X1;expand=cc1.*w;
w=sum(expand);
w1=cumsum(expand);
delta= sqrt(2*sum((w-y).^2)/nn);
disp(sprintf('> The RMS error is %g ', delta));

str='10';

gbutt([.75 .01 .15 .05], 'CONTINUE', 'uiresume');
tt2=text0([.15 .45 .75 .05],'Choose the number of terms in the expansion');
ee=edit1([.75 .45 .15  .05],str);
N=eval(ee);
delete(tt2),
rbutt([.75 .01 .15 .05], 'WAIT..', '');

for iter=1:N

	l1=plot(X,y,'-. r');axis(xy);hold on;
	set(l1,'LineWidth',2);
	plot(X,w1(iter,:));hold off;
	title(sprintf('Expansion up to %g terms',iter));	
	pause(1);
	
end

delta = sqrt(2*sum((w1(N,:)-y).^2)/nn);
xlabel(sprintf('The RMS error is %g ',delta))

q=0; %%%%%%%%%%%%%%%%%%%%%%%%%

bbutt([.60  .01  .15  .05],'More terms','uiresume; q=1;'); 
gbutt([.75  .01  .15  .05],'CONTINUE','close;q=0;'); 
uiwait;
N1=N;

while q==1 %%%%%%%%%%%%%%%%%%%%%%%%

N1= min([N1+5,N0]);

	l1=plot(X,y,'-. r');axis(xy);hold on;
	set(l1,'LineWidth',2);
	plot(X,w1(N1,:));hold off;
	title(sprintf('Expansion up to %g terms',N1));
	delta = sqrt(2*sum((w1(N1,:)-y).^2)/nn);
	xlabel(sprintf('The RMS error is %g ',delta))

bbutt([.60  .01  .15  .05],'More terms','uiresume; q=1;'); 
gbutt([.75  .01  .15  .05],'CONTINUE','close;q=0;'); 
uiwait;

if N0 <= N1
q=0;
end 

end

%%%%%%% show error term

	plot(X,w1(N1,:)-y),
	title(sprintf('This is the error term, RMS value = %g',delta)),
 	xlabel('Press a button to continue!');

   end

rbutt([.3   .01  .15  .05],'New f(X)','uiresume; q=2;');
rbutt([.45  .01  .15  .05],'REPEAT','close;q=3;'); 
bbutt([.6   .01  .15  .05],'BACK','close;q=4;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q=0;'); 
uiwait;

if q < 4 & q > 0

obutt,
		
elseif q==4 
		
special;return;
		
end

end

disp('> Type <legexp> to repeat this!');


