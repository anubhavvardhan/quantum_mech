
%> The file <well1.m> calculates the bound states in a number of different
%> quantum wells consisting of intervals with constant or linear potential.
%> The transfer matrix method is used, and an automatic search for the 
%> eigenvalues. The function files <trf1.m> and <trf2.m> are used to 
%> calculate the transfer matrices.
%>
%> © Goran Lindblad - gli@theophys.kth.se

clear, close;disp('> Welcome to <well1>!');

txt={' CALCULATION OF BOUND STATES IN A FAMILY OF QUANTUM WELLS' 
 ' ' 
 ' ' 
 ' The transfer matrix is calculated for a range of energies, and ' 
 ' the eigenvalues are found from the boundary conditions - i.e. ' 
 ' the continuity with the exponentially decaying solutions' 
 ' outside the interval. ' 
 ' ' 
 ' Choose the two parameters in the potential! ' };
 
p1=subplot(2,2,1);axis('off'),
q=get(p1,'Position');
l1=line([0,1],[1,1]);set(l1,'LineWidth',3),
l1=line([1,1],[0,1]);set(l1,'LineWidth',3),
l1=line([1,2],[0,0]);set(l1,'LineWidth',3),
l1=line([2,2],[0,1.5]);set(l1,'LineWidth',3),
l1=line([2,3],[1.5,1.5]);set(l1,'LineWidth',3),
gbutt([q(1)+0.3  q(2)  .05 .05],'#1','close;q0=1;'); 

%%%%%%%%%%%%	
p2=subplot(2,2,2);
q=get(p2,'Position');
axis('off');%axis('equal')
l1=line([0,1],[1,1]);set(l1,'LineWidth',3),
l1=line([1,2],[1, 0]);set(l1,'LineWidth',3),
l1=line([2,2],[0,1.5]);set(l1,'LineWidth',3),
l1=line([2,3],[1.5,1.5]);set(l1,'LineWidth',3),
gbutt([q(1)+.3  q(2)  .05 .05],'#2','close;q0=2;'); 

%%%%%%%%%%%%

p3=subplot(2,2,3);
q=get(p3,'Position');
% axis([0 1 0 1]);
axis('off');%axis('equal')

l1=line([0,.5],[1,1]);set(l1,'LineWidth',3),
l1=line([.5,.5],[1,0]);set(l1,'LineWidth',3),
l1=line([.5,1],[0,0]);set(l1,'LineWidth',3),
l1=line([1,1],[0,-.5]);set(l1,'LineWidth',3),
l1=line([1,1.5],[-.5,-.5]);set(l1,'LineWidth',3),
l1=line([1.5,1.5],[-.5,1]);set(l1,'LineWidth',3),
l1=line([1.5,2],[1,1]);set(l1,'LineWidth',3),
gbutt([q(1)+.3  q(2)  .05  .05],'#3','close;q0=3;'); 

%%%%%%%%%%%%
p4=subplot(2,2,4);
q=get(p4,'Position');
axis('off');%axis('equal')
l1=line([0,0.5],[1,1]);set(l1,'LineWidth',3),
l1=line([.5,.5],[1,0]);set(l1,'LineWidth',3),
l1=line([.5, 1],[0,0]);set(l1,'LineWidth',2),
l1=line([1,1],[0,1]);set(l1,'LineWidth',3),
l1=line([1,1.5],[1,1]);set(l1,'LineWidth',3),
l1=line([1.5,1.5],[1,-.5]);set(l1,'LineWidth',3),
l1=line([1.5,2],[-.5,-.5]);set(l1,'LineWidth',3),
l1=line([2,2],[-.5,1]);set(l1,'LineWidth',3),
l1=line([2,2.5],[1,1]);set(l1,'LineWidth',3),
gbutt([q(1)+.3  q(2) .05 .05],'#4','close;q0=4;'); 

suptitle('Bound states in quantum wells: choose your problem'),

uiwait; close,

subplot(2,1,1), axis('off'),
tt1 = text0([.15  .6 .75 .3], txt);

subplot(2,1,2),axis('off'),

if q0==1

l1=line([0,1],[1,1]);set(l1,'LineWidth',3),
l1=line([1,1],[0,1]);set(l1,'LineWidth',3),
l1=line([1,2],[0,0]);set(l1,'LineWidth',3),
l1=line([2,2],[0,1.5]);set(l1,'LineWidth',3),
l1=line([2,3],[1.5,1.5]);set(l1,'LineWidth',3);
t2= text(0.6,0.05,'- V_{0} < 0'); set(t2,'FontSize',18),
t2=text(1.65,1.5,'V_{1} > 0');set(t2,'FontSize',18),
t2=text(1.05,1.05,'X = 0');set(t2,'FontSize',18),
t2=text(2.05,1.05,'X = 1');set(t2,'FontSize',18),

elseif q0==2

l1=line([0,1],[1,1]);set(l1,'LineWidth',3),
l1=line([1,2],[1,0]);set(l1,'LineWidth',3),
l1=line([2,2],[0,1.5]);set(l1,'LineWidth',3),
l1=line([2,3],[1.5,1.5]);set(l1,'LineWidth',3);
t2= text(0.8,0.05,'- V_{0} < 0'); set(t2,'FontSize',18);
t2=text(1.65,1.5,'V_{1} > 0');set(t2,'FontSize',18);
t2=text(1.05,1.05,'X = 0');set(t2,'FontSize',18);
t2=text(2.05,1.05,'X = 1');set(t2,'FontSize',18);

elseif q0==3

l1=line([0,.5],[1,1]);set(l1,'LineWidth',3),
l1=line([.5,.5],[1,0]);set(l1,'LineWidth',3),
l1=line([.5,1],[0,0]);set(l1,'LineWidth',3),
l1=line([1,1],[0,-.5]);set(l1,'LineWidth',3),
l1=line([1,1.5],[-.5,-.5]);set(l1,'LineWidth',3),
l1=line([1.5,1.5],[-.5,1]);set(l1,'LineWidth',3),
l1=line([1.5,2],[1,1]);set(l1,'LineWidth',3),
t2= text(0.6,0.15,'- V_{0} + V_{1}'); set(t2,'FontSize',18),
t2=text(1.1,-.35,'- V_{0} < 0');set(t2,'FontSize',18),
t2=text(.4,1.1,'X = 0');set(t2,'FontSize',18),
t2=text(.9,1.1,'X = 1');set(t2,'FontSize',18),
t2=text(1.4,1.1,'X = 2');set(t2,'FontSize',18),


elseif q0==4

l1=line([0,0.5],[1,1]);set(l1,'LineWidth',3),
l1=line([.5,.5],[1,0]);set(l1,'LineWidth',3),
l1=line([.5, 1],[0,0]);set(l1,'LineWidth',2),
l1=line([1,1],[0,1]);set(l1,'LineWidth',3),
l1=line([1,1.5],[1,1]);set(l1,'LineWidth',3),
l1=line([1.5,1.5],[1,-.5]);set(l1,'LineWidth',3),
l1=line([1.5,2],[-.5,-.5]);set(l1,'LineWidth',3),
l1=line([2,2],[-.5,1]);set(l1,'LineWidth',3),
l1=line([2,2.5],[1,1]);set(l1,'LineWidth',3),
t2= text(0.58,0.15,'- V_{0} + V_{1}'); set(t2,'FontSize',18),
t2=text(1.55,-.35,'- V_{0}  < 0');set(t2,'FontSize',18),
t2=text(.4,1.1,'X = 0');set(t2,'FontSize',18),
t2=text(.9,1.1,'X = 1');set(t2,'FontSize',18),
t2=text(1.4,1.1,'X = 2');set(t2,'FontSize',18),
t2=text(1.9,1.1,'X = 3');set(t2,'FontSize',18),

end

str='[50,20]';
txt={'Set [V0,V1], both positive at right: [V0,V1]='
'Then multiply them using the sliders!'};
t2 = text0([.15  .5 .75 .1], txt);
gbutt([.75 .01  .15 .05], 'CONTINUE', 'uiresume');
VV=edit1([.75 .55 .15  .05],str);
VV=eval(VV);
delete(tt1), delete(t2),

str0='V0 = xx, V1 = yy. Change the values with the sliders!';

q1=2

while q1==2 %%%%%%%%%%%%%%%%%%%%%%

str=sprintf(strrep(str0,'xx','%g'),VV(1));  
str1=sprintf(strrep(str,'yy','%g'),VV(2));
ee=linspace(0,VV(1),700);
k1=sqrt(2*ee);

if q0==1

trans=trf1([0,1],-VV(1), - VV(1) + ee);
k2= sqrt(2*ee + 2*VV(2));

elseif q0==2

k2= sqrt(2*ee + 2*VV(2));
trans=trf2([0,1],[0,-VV(1)],- VV(1) + ee);

elseif q0==3

tran1=trf1([0,1],-VV(1)+VV(2), - VV(1) + ee);
tran2=trf1([1,2],-VV(1),- VV(1) + ee);
trans=mult4(tran1,tran2);
k2=k1;

elseif q0==4

tran1=trf1([0,1],-VV(1)+VV(2), - VV(1) + ee);
tran2=trf1([0,1],0,- VV(1) + ee);
trans=mult4(tran1,tran2);
tran3=trf1([0,1],-VV(1),- VV(1) + ee);
trans=mult4(trans,tran3);
k2=k1;

end

kk=[k2; ones(size(ee)); k1.*k2; k1];
wron=sum(trans.*kk);
wron1=crop(wron,100,-100);

zz=findzero(ee,wron);



if isempty(zz)

subplot(2,1,2), plot(- VV(1) + ee,wron),axis tight, 

str2='There are no bound states!';

t1 = text0([.15  .4 .75 .05], str2);
t2 = text0([.15  .5 .75 .05], str1);
hold off,

else

subplot(2,1,2), plot(-VV(1)+ee,wron1),axis tight, 
title('Eigenvalues are given by crossing of y=0 axis'),hold on,
plot(- VV(1) + ee, zeros(size(ee)), 'g -'),
plot(- VV(1) + zz, zeros(size(zz)), 'r o'),hold off;

nn=length(zz);
str=-VV(1)+zz';
disp('Eigenvalues = '), disp(str),

subplot(2,1,1), xy=[.5  nn+.5  -VV(1)  0];
t1= bar(-VV(1)+zz); axis(xy), title('Eigenvalues'),
t2 = text0([.15  .5 .75 .03], str1);

end

mult0=0; mult1=0;

g0 = slider0([.91  .15  .04  .75],[-2,2,mult0],'q1=2;uiresume;');
g1 = slider0([.96  .15  .04  .75],[-2,2,mult1],'q1=3;uiresume;');

rbutt([.45  .01  .15  .05],'REPEAT','uiresume, q1=1;'),
bbutt([.6   .01  .15  .05],'BACK','uiresume,q1=0;'); 
bbutt([.75  .01  .15  .05],'QUIT','close;q1=-1;'); 
uiwait; 

if q1==2
mult0=get(g0,'Value');VV(1)=VV(1)*exp(mult0);delete(t1,t2);
elseif q1==3
mult1=get(g1,'Value');VV(2)=VV(2)*exp(mult1); q1=2;
delete(t1,t2);
end 

end  %%%%%%%%%%%% end of while q1 = 2 

if q1==1
well1; return;
elseif q1==0
boundst, return;
end
close;
disp('Type <well1> to do this again!');


