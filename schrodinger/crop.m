function f = crop(y, ymax, ymin);
%> The function file <crop.m> crops a real vector or matrix: 
%> it replaces values y > ymax by ymax, values < ymin by ymin. 
%> Call: crop(y, ymax, ymin);
%>
%> © Goran Lindblad - gli@theophys.kth.se

% First make sure ymin < ymax

if ymax < ymin

z=[ymin,ymax];

ymax=z(1); ymin=z(2);

disp('> Interchanging ymax and ymin!');

elseif ymax == ymin

disp('> You have chosen ymax = ymin, a constant function is returned!');
 
end

% Now do the cropping 

wmax= 0.5*(sign(y-ymax) + 1); % 1 if y > ymax, 0 otherwise.
wmin= 0.5*(-sign(y-ymin) + 1); % 1 if y < ymin, 0 otherwise.
f=y.*(1-wmax).*(1-wmin) + ymax*wmax + ymin*wmin;

