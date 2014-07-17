%> <well3.m> solves for the eigenvalues of a quadratic and a circular 
%> infinitely deep well in 2D, of the same area, and compares the spectra.
%> Variable parameters, A = area, E = the maximal energy
%> The search for eigenvalues i automatic, using <findzero.m>
%>
%> © Goran Lindblad - gli@theophys.kth.se

%% GL 990428

clear, close; disp('> Welcome to <well3>');

txt={' EIGENVALUES OF AN INFINITE 2D WELL ' 
' ' 
' The script <well3.m> solves for approximate eigenvalues for a 2D ' 
' quantum well of quadratic or circular shape of given area A and  '
' infinitely deep.' 
''
' It solves for the eigenvalues in a range [0,E] from the boundary '
' condition: the solutions regular at the origin are set to zero on '
' the boundary. '
''
' The quadratic case is solved by first solving the 1D problem, and '
' then adding the resulting spectra to the spectrum of the 2D case.'
' In the circular well the radial functions are Bessel functions of'
' type J_m.  The numerical search for eigenvalues proceeds for '
' increasing m = 0,1,... until no more eigenvalues are found. '
' The script <besjz.m> is used to find them.'
''
' You can choose the values of A and E.'
''
' '
''
' Insert new values of [A,E] in the box >> '
' or press the continue button!'};

tt1=text0([.12 .3 .78 .6], txt);
gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume');
str='[ 20 , 20 ]';
ee=edit1([.75 .3 .15  .05],str);
AE=eval(ee); A=AE(1); Emax=AE(2);

delete(tt1);

xdim= sqrt(A); rdim=sqrt(A/pi);

nmax=sqrt(2*Emax)*xdim/pi;

% nn=100; % number of eigenvalues
nx=1:nmax;

%% Energies for 1D infinite well
kk= nx*pi/xdim; 
e1=kk.^2/2; 

%% Find the energies for the  2D infinite well
[x,y]=meshgrid(e1,e1);
e2=x+y;
e2=reshape(e2,length(e1)^2,1);
e2=sort(e2);

test=e2-Emax;
test=crop(test,0,-Emax);
ind=find(test);
e2=e2(ind);


bar(e2);axis([.5, length(e2)+.5, 0, Emax]),
title('Bar diagram of the lowest states, quadratic case');
ylabel('Energy eigenvalue'),
xlabel('Ordinal of the state');
cbutt(.75, .01);

l2=length(e2);
[x,y]=stairs(e2);
x=[0;0;x-1;l2;l2];y=[0;0;y;e2(l2);e2(l2)];
p1=plot(y,x); axis([0, e2(l2)+.5, 0, l2+.5]),
title('Distribution function for eigenvalues, quadratic case'),
xlabel('Energy E');
ylabel('Number of eigenstates with energy <= E');
set(p1,'LineWidth',3),

cbutt(.75, .01),

q1=1; m=0;
eee=[];
rbutt([.75  .01  .15  .05],'WAIT...',''); drawnow;

while q1==1

nmax= ceil((.5 + sqrt(2*Emax)/pi )*rdim);

zr=besjz(m,nmax);

er=zr.^2/2/rdim/rdim;

test=er-Emax;
test=crop(test,0,-Emax);
ind=find(test);

if isempty(ind)

q1=0;

else
er=er(ind);

if m > 0  %% double degeneracy!

er=[er;er];
er=sort(er);

end 

eee=mad(eee,er');

end


m=m+1;

% pause(1)
% cbutt(.75, .01)


end 


le=size(eee);

bar([0:le(1)-1],eee,'b');
axis([-.5 , le(1)-.5, 0, Emax]);
title('Circular case: energy spectra for  m = 0, 1, ..');
xlabel('Angular momentum m');
cbutt(.75, .01)


eee=reshape(eee,le(1)*le(2),1);
eee=sort(eee);

ind = find(eee);
er=eee(ind);
bar(er); axis([.5, length(er)+.5, 0 ,Emax]);
title('Bar diagram of energy spectrum,circular case');
cbutt(.75,.01)

lr=length(er);
[xr,yr]=stairs(er);
xr=[0;0;xr-1;lr;lr];yr=[0;0;yr;er(lr);er(lr)];
p1=plot(yr,xr); axis([0, er(lr)+.5, 0, lr+.5]),
title('Distribution function for eigenvalues, circular case'),
xlabel('Energy E');
ylabel('Number of eigenstates with energy <= E');
set(p1,'LineWidth',3),

cbutt(.75,.01),
hold on,
p1=plot(y,x,'r'); set(p1,'LineWidth',3),
title('Distribution function for eigenvalues, circular and quadratic cases'),
hold off,
rbutt([.3   .01   .15  .05],'REPEAT','close; q1=1;'),
bbutt([.45  .01   .15  .05],'BACK','close;q1=2;') 
bbutt([.6   .01   .15  .05],'MAIN MENU','close;q1=3;'), 
bbutt([.75  .01   .15  .05],'QUIT','close;q1=4;'), 
uiwait; 

if q1==1
well3;
elseif q1==2
close; inbox; return;
elseif q1==3	
close; start; return;
else
close;
disp('> Type <well3> to do this again!');
end 

%%%%%%%%%%
