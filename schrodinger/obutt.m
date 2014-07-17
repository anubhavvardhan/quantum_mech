function h = obutt;
%> The file <obutt.m> gives a grey non-functional button
%> at a predetermined position at the bottom of the figure.
%> Call: obutt
uicontrol('Style','pushbutton','Units','normalized',...
          'Position',[.15 .01 .75 .05],'BackGroundColor',[.8 .8 .8])
%>
%> © Goran Lindblad - gli@theophys.kth.se
