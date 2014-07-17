function DY=sde(X,Y);
%> The file <sde.m> defines a Schrodinger type 1D time-independent DE
%> with a potential V defined by a string <potential> 
%> which is a global variable with space variable X. 
%> For instance: <potential = '15*pot1(X)'>  points to the 
%> function file <pot1.m>.  The value of the energy is given
%> by the global variable <energy>, e.g. : <energy = 12>. 
%> The normalization is hbar = mass = 1.
%> Call: [X,Y]=ode45('sde',[0,1],[y;Dy]);
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961125

global potential energy

V=eval(potential); E=energy;

DY1=Y(2); DY2=2*(V - E).*Y(1);

DY=[DY1;DY2];

 
%%% © Göran Lindblad 1996
