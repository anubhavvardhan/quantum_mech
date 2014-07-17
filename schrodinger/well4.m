%> <well4.m> solves for the eigenvalues of a rectangular, infinitely 
%> deep well in 2D, of variable area and ratio. It calculates the
%> spectrum and the spectrum of energy differences.
%> Variable parameters:
%> A = area, R = ratio of rectangle, E = the maximal energy
%> The search for eigenvalues i automatic, using <findzero.m>
%> © Goran Lindblad - gli@theophys.kth.se

clear, close; disp('> Welcome to <well4>');

q1=1;

% starting values for A, R, E
xyz='[50, 3, 30]';

txt={' EIGENVALUES OF AN INFINITE 2D WELL ' 
' ' 
' The script <well4.m> solves for the eigenvalues for a 2D infinitely' 
' deep quantum well of rectangular shape. The variable parameters are'
' the area A, the aspect ratio R, and maximal energy E. '
''
' It solves for the eigenvalues from the boundary condition that the'
' solutions are zero on the boundary. '
''
' The spectrum is found by first solving the two 1D problems, and '
' then adding the resulting spectra to the spectrum of the 2D case.'
' The distributions of the spectrum and of the energy differences'
' are calculated. Note that the spectrum of energy differences is'
' sensitive to the rationality/irrationality of the aspect ration!'
'' 
' You can change the values of A, R, and  E. Depending on your '
' choice the spectrum may turn out to be empty!' 
''
''
''
' Now insert new values of [A, R, E] in the box >> '
' or press the continue button!'};

tt1=text0([.12 .3 .78 .6], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
xyz=edit1([.75 .3 .15  .05],xyz);
ee=eval(xyz); A=ee(1); R=ee(2);  Emax=ee(3);
delete(tt1);

while q1==1

xdim= sqrt(A); %% side of square of area A
rr=sqrt(R); %% scaling factor
%% max no eigenvalues for 1D system 
nmax=sqrt(2*Emax)*xdim/pi;

%% Now scale sides, to get a rectangle of given A,R from
%% square of same area

x1=xdim*rr; x2=xdim/rr; 
n1=ceil(rr*nmax); n2=ceil(nmax/rr); %% scaled no of eigv
nn1=1:n1; nn2=1:n2;

%% Now solve for the two 1D problems
%% solutions sin(kx) with k xdim = 0, k > 0

kk1= nn1*pi/x1; 
e1=kk1.^2/2; % energies for 1D infinite well
%% cut out values larger than Emax
test=crop(e1-Emax,0,-Emax); 
ind=find(test);
e1=e1(ind);

%% you may end up with no eigenvalues at all
if isempty(ind);
disp('> The spectrum is empty!')
q1=0;

end

%% Now for the other 1D system

if q1==1

kk2= nn2*pi/x2;
e2=kk2.^2/2; % energies for 1D infinite well
test=crop(e2-Emax,0,-Emax);
ind=find(test);
e2=e2(ind);

if isempty(ind);
disp('> The spectrum is empty!')
q1=0;

end

end 

%% Now sum the two spectra, order them, and cut out the 
%% values larger than Emax

if q1==1

[x,y]=meshgrid(e1,e2);
ee2=x+y;
ee2=reshape(ee2,length(e1)*length(e2),1);
ee2=sort(ee2);

test=crop(ee2-Emax,0,-Emax);
ind=find(test);
ee2=ee2(ind);
ll2=length(ee2);
nn2=[1:ll2];

subplot(2,1,1)
bar(e1), axis tight,
title('Spectra of the two 1D systems');
subplot(2,1,2)
bar(e2), axis tight,

cbutt(.75 ,.01);

txt=sprintf('The values are A =%g, R = %g, E = %g',A,R,Emax);

subplot(1,1,1)
plot(nn2,ee2);axis tight,
tt1=text0([.15 .85 .5 .05], txt);
title('Spectrum for the rectangular infinite well');
cbutt(.75 ,.01);

[x2,y2]=stairs(ee2);
x2=[0;0;x2-1;ll2;ll2];y2=[0;0;y2;ee2(ll2);ee2(ll2)];
p1=plot(y2,x2); 
axis([ee2(1)-.5, ee2(ll2)+.5, 0, ll2+.5]),
title('Distribution function for energy eigenvalues');
cbutt(.75 ,.01);

de2=diff(ee2);
de2=sort(de2);
stairs(de2); axis tight,
title('Spectrum of energy differences');
cbutt(.75 ,.01);
delete(tt1),

end




if q1==0

plot(linspace(0,1), 0*linspace(0,1)); axis off,
title('The spectrum is empty!!');
q1=1;

end 

txt={' Now insert new values of [A,R, E] in the '
' box or press the continue button!'};
tt1=text0([.15 .8 .6 .1], txt);
gbutt([.75  .01  .15  .05],'PASS','q1=0;uiresume'); 
xyz=edit1([.6 .8 .15  .05],xyz);
ee=eval(xyz); A=ee(1); R=ee(2);  Emax=ee(3);
delete(tt1);


end
str='> Type <well4> to do this again!';
rbutt([.3   .01   .15  .05],'REPEAT','uiresume;well4;'); 
bbutt([.45  .01   .15  .05],'BACK','inbox;'); 
bbutt([.6   .01   .15  .05],'MENU','start;'); 
bbutt([.75  .01   .15  .05],'QUIT','close; disp(str)'); 
uiwait;

return

