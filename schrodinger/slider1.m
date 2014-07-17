function f = slider1(xy,rr)
%> call function slider1(xy,rr)
%> xy = position 4-vector [x1 y1 dx dy]
%> rr = value  3-vector [min, max, value]
q= uicontrol('Style','slider','Units','normalized',...
'Position',xy,'BackGroundColor',[.4 .2 .4],'ForeGroundColor','w',...
'Min',rr(1),'Max',rr(2),'Value',rr(3),...
'Callback','uiresume');
uiwait;
h=get(q,'Value');
delete(q);
f=h;
%>
%> © Goran Lindblad - gli@theophys.kth.se
