
%> The file <choice.m> allows you to choose the potential in the 
%> Schrodinger equation from the list in the function file <pot1.m>.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125
clear, close,

global PP pot1str 

q1=1;

while q1==1

q0=1; string='uiresume;q0=P;';
x=linspace(0,1,200); zx=zeros(size(x)); ux=ones(size(x));

	for PP=1:11

		B=[pot1(x);zx];
		p=subplot(4,3,PP);  	%fill(x,B,'r')
		q=get(p,'Position');
		q=[q(1)+0.0,q(2)+0.1,0.05,0.05];
		P1=strrep(string,'X','-X'); %
		string0=sprintf(strrep(string,'P','%g'),PP);
		hold on;
		gbutt(q,PP,string0); 
		fig=plot(x,B,'r');
		axis('off');

		plot(x,zx,'k');
		hold off
	end
	
	for PP=12

	p=subplot(4,3,PP);  	%fill(x,B,'r')
	q=get(p,'Position');
	q=[q(1)+0.0  q(2)+0.1 0.05 0.05];
	txt={	' User defined by'
			' string input'};
	t1=text0([q(1) q(2)-.08 .2 .07],txt);
	string0=sprintf(strrep(string,'P','%g'),PP);
	hold on;
	gbutt(q,PP,string0); 
 	axis('off');
	
	end

suptitle('Choose the potential by pressing the button!');
uiwait;
disp(sprintf('> You have chosen # %g',q0));
PP=q0;

	if PP==12
	
	txt={'Choose the potential as a string in the window!'};
	pot1str='-1 + 2*abs(X-.5)';
	t1=text0([.15 .4 .75 .3 ],txt);
	gbutt([.75 .01 .15 .05], 'CONTINUE','uiresume'),
	pot1str=edit1([.15 .5 .75 .05],pot1str);
	delete(t1);
	
	end

G=200; % no of points in display of potential and wave functions
X=linspace(0,1,G);
V=pot1(X);% V=eval(pot1str);

%% Display of potential:
close,
l1=plot(X,V); axis('equal'), axis('tight'), set(l1,'LineWidth',3),
title('This is your potential!'),
if PP==5
txt={' With a multiplier - 32*n*(n+1), n = 1,2,.. this potential is exactly'
	' solvable!  The eigenvalues are then - 32*k^2 , k= 0,1,..,n.'
	' It is also reflectionless, the tranmission coefficient = 1!'};
t1=text0([.15 .6 .75 .1 ],txt);
end

rbutt([.6  .01  .15  .05],'CHANGE','close, q1=1;'); 
gbutt([.75  .01  .15  .05],'CONTINUE','close,q1=0;');
uiwait;
end 

