function f = slider0(xy,rr,cc)
%> call function slider0(xy,rr,cc)
%> xy = position 4-vector [x1 y1 dx dy]
%> rr = value  3-vector [min, max, value]
%> cc = callback string eg 'uiresume'
q = uicontrol('Style','slider','Units','normalized',...
'Position',xy,'BackGroundColor',[.4  .2 .4],'ForeGroundColor','w',...
'Min',rr(1),'Max',rr(2),'Value',rr(3),...
'Callback',cc);
f=q;

%>
%> © Goran Lindblad - gli@theophys.kth.se
