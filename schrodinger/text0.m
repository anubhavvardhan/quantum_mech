function f = text0(xy,ww)
%> Call f=text0(xy,'string');
%> where x is a 4 - vector [x0 y0 dx dy];
%> ww is a string in cell array .
f=uicontrol('Style','text','Units','normalized','FontSize',12,...
'String',ww, 'HorizontalAlignment','left','Position',xy,...
'BackgroundColor',[.93 .93 .93]);
