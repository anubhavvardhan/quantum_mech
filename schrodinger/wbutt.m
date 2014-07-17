function h = wbutt(x,y)
%> This is a standard red 'wait' button of given 
%> length = .15 and height = .05.
%> Call:  wbutt(x,y)
%> Input: x, y = position of lower left corner
%> No callback.

uicontrol('Style','pushbutton','Units','normalized',...
          'Position',[x,y,.15, .05],'String','WAIT...',...
          'BackGroundColor',[.85 .0 .0],'ForeGroundColor','w',...
          'Fontsize',12); 

%>
%> © Goran Lindblad - gli@theophys.kth.se
