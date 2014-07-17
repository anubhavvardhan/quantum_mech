function f = mult4(x,y)
%> The function file <mult4>  multiplies two 2x2 matrices represented
%> as column 4-vectors, where each matrix element is a row vector, 
%> of same dimension for x and y.
%> Call: mult4(x,y)
%> Input: x, y are 4 x N matrices,
%> Output: a matrix of same dimension

f=[ x(1,:).*y(1,:) + x(2,:).*y(3,:);
    x(1,:).*y(2,:) + x(2,:).*y(4,:);
    x(3,:).*y(1,:) + x(4,:).*y(3,:);
    x(3,:).*y(2,:) + x(4,:).*y(4,:)]; 
    
%>
%> © Goran Lindblad - gli@theophys.kth.se
