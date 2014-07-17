function h = cbutt(x,y)
%> This is a standard green 'continue' button of given 
%> length = .15 and height = .05.
%> Call:  cbutt(x,y)
%> Input: x, y = position of lower left corner
%> The callback is 'uiresume'.

uicontrol('Style','pushbutton','Units','normalized',...
          'Position',[x,y,.15, .05],'String','CONTINUE',...
          'BackGroundColor',[.0 .7 .0],'ForeGroundColor','w',...
          'Fontsize',12,'Callback','uiresume'); uiwait;

%>
%> © Goran Lindblad - gli@theophys.kth.se
