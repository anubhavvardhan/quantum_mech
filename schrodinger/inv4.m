function f = inv4(x)
%> This file produces the inverse of one 2x2 matrix represented as a 
%> column 4-vector. Each matrix element is a row N-vector.
%> Call: inv4(x)
%> Input: x = 4 x N matrix.
%> Output: a matrix of same dimension.
%>
%> © Goran Lindblad - gli@theophys.kth.se
 
det = [1,1,1,1]'*(x(1,:).*x(4,:)-x(2,:).*x(3,:));

f=[ x(4,:);
   -x(2,:);
   -x(3,:);
    x(1,:)]./det;
    

