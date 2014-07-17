function f=glue(y,z);
%> <glue.m> glues jumps in given functions producing  'continuous' 
%> functions as output. This assumes that the jumps in function values
%> are multiples of a positive number z which is given as input. 
%> Examples are jumps of value \pi in phase shifts etc.
%> Call glue(y,z); 
%> y a M x N matrix, each row is a function  y=y(x) which is to
%> be glued,
%> z > 0 is the smallest common factor in the jumps.
%>
%> © Goran Lindblad - gli@theophys.kth.se

% GL 961123

x=size(y);N=x(1);

yy=y/z; y1=diff(yy'); 

y1=[zeros(1,N);round(y1)]; y1=cumsum(y1);yy=yy-y1';

f=z*yy;
