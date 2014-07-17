function f=mad(A,x);
%> The function file <mad.m> adds a row vector to an existing matrix
%> Call: mad(a,x)
%> Input: a is a matrix, x a row vector.
%> Output: matrix with x as a new last row to a, 
%> expanding the dimensions as needed.
%>
%> © Goran Lindblad - gli@theophys.kth.se

%  GL 961125

if isempty(A)

f=x;

elseif isempty(x)

f=A;

else
sA=size(A);
sx=size(x);
if sx(1) > 1
disp('The vector x should be a row vector!');
return
end 
if sx(1)==0
return
end
if sA(2)==sx(2)
A=[A;x];
end
if sA(2) > sx(2)
x=[x,zeros(1,sA(2)-sx(2))];
A=[A;x];
end
if sA(2) < sx(2)
A=[A,zeros(sA(1),sx(2)-sA(2))];
A=[A;x];
end

f=A;

end

%%% © Göran Lindblad 1996
