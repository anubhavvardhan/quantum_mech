function f=edit1(xy,string);
%> The file <edit1> allows us to edit a string
%> Call q = edit1(xy, string) where
%> xy = 4-vector giving the position and size of window
%> string = string to be edited. 
%> If you want a number, say q=eval(q) 

t0 = uicontrol('Style','edit','Units','normalized', ...
	'BackgroundColor',[.9 .9 .9], ...
	'Position',xy, 'String',string, ...
	'Callback','uiresume');
uiwait; 

f=get(t0,'String');delete(t0);

%>
%> © Goran Lindblad - gli@theophys.kth.se
