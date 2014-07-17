function WPACKET(Command)
%	PROGRAM WPACKET.M, GUI-INTERFACE for MATLAB 5 
%
%	solves the 1-dimensional Schrödinger equation for various
%	potentials and parameters. The routine calls the following files:
%	
%	GSSPCKT.M	(generates a Gaussian wave-packet)
%	PTNTLS.M	(a collection af various energy potential profiles)
%	SCHROEQ.DLL   (core routine performing the time evolution of wave amplitude)
%	

global Psi x xWidth pNum runFlag stopFlag Comment V xo Dx
 
if nargin < 1,
    Command = 'Initialize';
	figure(1)
	clf
    set(1, 'Name', ['Schrödinger wave packets in 1 dimension',...
	       '  -  R.Nobili - Padova University - ITALY']);         
	set(1, 'Userdata', zeros(size(1:30)));
 end
          
if ~strcmp(Command, 'Initialize'),
   handles=get(1, 'Userdata'); % retrieve data saved at the end of
   % the 'Initialize' call; used by all other calls:      
   H_AX=handles(1); % save data in the Userdata buffer of the window
   H_LV=handles(2);
   H_LP=handles(3);
   H_LR=handles(4);
   H_LI=handles(5);
   H_EkTxt=handles(6);
   H_LEav=handles(7);	
   H_TXT=handles(8);
   % H_TXT0 =handles(9);       
   H_POP=handles(10);
   % H_FIG=handles(11);
   H_TEav=handles(12);
	H_EavTxt=handles(13);
  	relMass=handles(14);
   % potFlag=handles(15);
	Ko=handles(16);
   Dx=handles(17);
   H_SL=handles(18);     
   H_ED=handles(19);
   Eav=handles(20);
    Ek=handles(21);
   lambda=handles(22);
   H_RUN=handles(23);
   H_STOP=handles(24);
   H_REST=handles(25);
   H_LamTxt=handles(26);
   H_DxTxt=handles(27);
   H_WED = handles(28);
   H_WSL = handles(29);
end

if strcmp(Command, 'Initialize'),
       
    % Default initialization:
    % handles=zeros(1,36);
   
	relMass=1; % default value in GSSPCKT.M
	xWidth=40;  % (Angstrom) default value in GSSPCKT.M
    pNum=1001;  % default value GSSPCKT.M
    Ko=0;	% initially the particle is at rest

	N=(pNum-1)/2;
	x=xWidth*(-N:N)/pNum;
	x=x(:);

	[V, xo, Dx, Comment] = PTNTLS(x,0); % load the potential profile list.
    
	[Psi, x, Ek, lambda] =  GSSPCKT(xo, Dx, Ko, relMass, xWidth, pNum);
    % generates the basic Gaussian packet as initial waveform configuration 
    
	Eav = Ek+sum(V.*(abs(Psi).^2)); % mean energy of the initial state
   
    potFlag=0; % zero if no potential profile is loaded
    runFlag=0; % zero if the waveform evolution is stopped

% ---------------- USER INTERFACE CONTROLS ------------------------------
		H_FIG=1;    	     	
    
  		H_RUN=uicontrol(1, 'style', 'pushbutton',...
         'Units', 'normalized',...
         'Position', [0.82,0.6,0.08,0.04],... % run button
			'string', 'Run','callback', 'WPACKET(''Run'');');       
 
		H_STOP=uicontrol(1, 'style', 'Pushbutton',...
         'Units', 'normalized',...
         'Position', [0.91,0.6,0.08,0.04],... % stop button
         'string', 'Stop',...
         'Enable', 'off', ...
         'callback','WPACKET(''Stop'');'); 

	    H_REST = uicontrol(1, 'Style', 'PushButton',... % restart button
         'Units', 'normalized',...
         'Position', [0.82,0.55,0.08,0.04],... 
         'string', 'Restart',...
         'Callback','WPACKET(''Rest'');',...
         'Enable', 'off');  
      
        txtEdPos=[0.0,0.0,0.3,0.04];
          
        H_TXT = uicontrol(1,...
                'Style', 'text',...
                'Units', 'normalized',...
                'Position', txtEdPos,...
                'FontSize', 10,... 	
                'String', ['Mean wave-number [1/',197,'] = '],...
                'Backgroundcolor', [0.8,0.8,0.8],...
                'HorizontalAlignment', 'right');    

        Ext=get(H_TXT, 'Extent');
       
        edPos=[0.01+Ext(3),0.0, 0.05,0.04];	
	
        H_ED = uicontrol(1,... % Edit window for wave-number width control
               'Style', 'edit',...
               'Units', 'normalized',...
               'Position', edPos,...
               'String', '0',...
               'HorizontalAlignment', 'center',...
               'Callback','WPACKET(''WNumberEditChanged'');',...  
               'Backgroundcolor', [1,1,1]);
          
		slPos = [Ext(3)+0.06, 0.0, 0.12, 0.04];	
        
        H_SL=uicontrol(1,'Style', 'slide',... % Slider for wave-number control
               'Min', 0,...
               'Max', 10,...
               'Value', 0,...   
               'Callback', 'WPACKET(''WNumberSliderMoved'');',... 
               'Units', 'normalized',...
               'Position', slPos);  
          
      txtWPos=[0.51,0.0,0.3,0.04];
    
      %H_WTXT = 
        uicontrol(1,...
			  'Style', 'text',...
             'Units', 'normalized',...
             'Position', txtWPos,...
	      	 'FontSize', 10,... 	
             'String', ['Gaussian semi-width [',197,'] = '],...
	      	 'Backgroundcolor', [0.8,0.8,0.8],...
             'HorizontalAlignment', 'right');    

       WExt=get(H_TXT, 'Extent');
       
       edWPos=[0.5+Ext(3),0.0, 0.08,0.04];	
	
       H_WED =uicontrol(1,...	% Edit window for Gaussian width control
             'Style', 'edit',...
             'Units', 'normalized',...
             'Position', edWPos,...
             'String', '0.9864',...
             'HorizontalAlignment', 'center',...
             'Callback','WPACKET(''GWidthEditChanged'');',...  
             'Backgroundcolor', [1,1,1]);
          
		slWPos=[0.51+WExt(3)+0.07, 0.0, 0.12, 0.04];	
        
        H_WSL=uicontrol(1,'Style', 'slide',... 	% Slider for Gaussian width control
             'Min', 0.2,...
             'Max', 10,...
             'Value', Dx,...   
             'Callback', 'WPACKET(''GWidthSliderMoved'');',... 
             'Units', 'normalized',...
             'Position', slWPos);  
         
		uicontrol(1, 'style', 'Pushbutton',...
            'Units', 'normalized',...
            'Position', [0.69,0.049,0.30,0.04],... % stop button
            'string', 'Reset Gaussian profile to default',...
            'Enable', 'on', ...
            'callback','WPACKET(''Res2Def'');'); 
 
% -------------- DATA PLOTTING ----------------------------
		uicontrol(1, 'style', 'text',...
         'string','Potential type', ...
         'Units', 'normalized',...
         'Position', [0.8, 0.94, 0.18, 0.03],...
         'Background', [0.8, 0.8, 0.8]);
      
        % LIST OF POTENTIAL PROFILES
	    potString = ['Zero potential well|High narrow barrier|Small barrier|',...
         'Medium barrier|Wide barrier|Large barrier|',...
         'Deep well|Small well|'...   
         'Medium well|Wide well|Large well|',...
         'Positive linear|Negative linear|',...
         'Positive periodic|Negative periodic|',...
         'Harmonic oscillator|Harm. narrow packet|',...
         'Harm. wide packet|Quadric symmetric|',...
         'Quadric tilted|Quadric strongly tilted|Special',...
         'Gaussian random'];
   
	    H_POP = uicontrol(1, 'style', 'pop','String',potString,...
        	'Callback', 'WPACKET(''Choose_pot'')',...
            'Value', 1,...
            'Units', 'normalized',...
            'Position', [0.8, 0.9, 0.2, 0.03],...
            'Background', 'w',...
            'Enable', 'on');   

        H_AX = axes('Units', 'normalized',...
            'Position', [0.075, 0.15, 0.725, 0.78], ...
            'XLim',[-xWidth/2,xWidth/2], 'XLimMode','manual',...
            'YLim',[-200,400],'YLimMode','manual', 'Box', 'on');   

        % MAIN WINDOW LEGENDS:
        text(0.75, 0.15, '- Probability density',...
            'Unit', 'Normalized',...
            'Color','k', 'FontName', 'Arial',...
            'FontSize', 8);
      
        text(0.75, 0.1, '- Real part',...
            'Unit', 'Normalized',...
            'Color',[1,0.5,1], 'FontName', 'Arial',...
            'FontSize', 8);

        text(0.75, 0.05, '- Imaginary part',...
            'Unit', 'Normalized',...
            'Color',[0, 0.8, 0.8], 'FontName', 'Arial',...
            'FontSize', 8);     
      
        xlabel('X-axis position (Angstrom)');
	
        % TITLE
        title('Potential Energy Profile in eV (thick blue line) & Wave Forms ');
	
        % BLUE LINE REPRESENTING THE POTENTIAL PROFILE
        H_LV =line('XData', x, 'Ydata', V, 'Color',...
            'b', 'LineWidth', 2, 'EraseMode', 'xor');

        % TEXT INDICATING THE MEAN ENERGY LEVEL
        H_TEav=text(-19, Eav+16,'<E>',...
            'EraseMode', 'background','Color', 'r',...
            'FontName', 'Arial','FontSize', 9);
        
        % PROBABILITY DENSITY PROFILE
        H_LP=line('XData', x,'YData',zeros(size(x)), 'Color',...
                'k','EraseMode', 'xor');
            
        % IMAGINARY PART OF WAVE PROFILE
        H_LR=line('XData', x, 'YData', zeros(size(x)),...
		'Color', 'c','EraseMode', 'xor');
    
        % REAL PART OF WAVE PROFILE
        H_LI=line('XData', x, 'YData', zeros(size(x)),...
            'Color', [1, 0.75, 1],'EraseMode', 'xor');
        
        % MEAN ENERGY LEVEL
        H_LEav=line('XData', [-20,20], 'YData', [Eav,Eav],...
            'Color', 'r','EraseMode', 'xor');
   
        % TEXT INDICATING THE ELAPSED TIME
        H_TXT0=text(0.02,0.05,'Time [fsec] = ','Units', 'normalized',...
      		'Horizontalalignment', 'left',...
      		'EraseMode', 'xor',...
      		'FontName', 'Arial',...
      		'FontSize', 8);

        Extent=get(H_TXT0, 'Extent');
        H_TXT =	text(0.02+Extent(3),0.05,'0.00','EraseMode', 'xor',...
      		'Horizontalalignment', 'left',...
   			'Units', 'normalized',...
      		'FontName', 'Arial',...
            'FontSize', 8);
         
        
        % DATA WINDOW % H_ASMALL = 
          axes('Units', 'normalized',...
            'Position', [0.805, 0.15, 0.195, 0.2], ...
            'Box', 'on',...
            'XLim',[0,1], ...
            'YLim',[0,1],...
            'FontName', 'Arial',...
            'FontSize', 9,...
            'YTickMode', 'manual',...
            'XTickMode', 'manual');
      
            title('Parameters at t=0');
    
        % DATA WINDOWS INSCRIPTIONS>
        H_EavTxt= text(0.02, 0.85, ['<E> = ',num2str(Eav),' eV'], ...
            'FontName', 'Arial','FontSize', 9); 
        H_EkTxt= text(0.02, 0.65, ['<E_{kin}> = ',num2str(Ek),' eV'], ...
            'FontName', 'Arial','FontSize', 9);
        H_LamTxt= text(0.02, 0.45,...
            ['<\lambda> = ', num2str(lambda),' ', 197], ...
            'FontName', 'Arial','FontSize', 9);
        H_DxTxt= text(0.02, 0.25, ['<\Delta x> = ', ...
            num2str(Dx),' ', 197], 'FontName', 'Arial','FontSize', 9);
       
%--------------- END DATA PLOTTING --------------------------------
%--------- SAVING DATA IN THE FIGURE-WINDOW USERDATA STACK handles -----------------	
  	
        handles=[H_AX, H_LV, H_LP, H_LR, H_LI,...
            H_EkTxt, H_LEav, H_TXT, H_TXT0,...
            H_POP, H_FIG, H_TEav, H_EavTxt, relMass,...
            potFlag, Ko, Dx, H_SL, H_ED, Eav, Ek, lambda,...
            H_RUN, H_STOP, H_REST, H_LamTxt, H_DxTxt, H_WED, H_WSL];
    
        set(1, 'Userdata', handles); % save data in the figure stack
 
        WPACKET('again');
   
 elseif strcmp(Command, 'WNumberSliderMoved'),
   
      if runFlag > 0 ,
 	     	set(H_SL,'value',Ko);
         return
       end
       Ko = get(H_SL,'value');  
       set(H_ED, 'String', num2str(Ko));
       handles(16)=Ko;  
       set(1,'Userdata', handles); 
	   WPACKET('again');
       
     
elseif strcmp(Command, 'WNumberEditChanged'),          
        if runFlag >0, ... 
		    Ko=handles(16);
 		    set(H_ED,'String',num2str(Ko));
            return;
        end;
        vStr = get(H_ED, 'String');
        Ko = str2double(vStr);
        
        if(Ko<0),
           Ko=0;
        end
        
        if(Ko>10),
           Ko=10;
        end
        
        set(H_ED, 'String', num2str(Ko));  
        set(H_SL, 'Value', Ko);
        handles(16)=Ko;
        set(1,'Userdata', handles); 
        WPACKET('again');

elseif strcmp(Command, 'GWidthSliderMoved'),
   
        if runFlag > 0,
 	     	set(H_WSL,'value',Dx);
        return
        end
        Dx = get(H_WSL,'value');  
        set(H_WED, 'String', num2str(Dx));
        handles(17)=Dx;
    
        set(1,'Userdata', handles); 
	    WPACKET('again');
       
elseif strcmp(Command, 'GWidthEditChanged'),
          
        if runFlag >0, ... 
		     Dx=handles(17);
 		     set(H_WED,'String',num2str(Dx));
             return;
        end;
        vStr = get(H_WED, 'String');
        Dx = str2double(vStr);
        
        if(Dx<0.2),
           Dx=0.2;
        end
        
        if(Dx>10),
           Dx=10;
        end
        
        set(H_WED, 'String', num2str(Dx));  
        set(H_WSL, 'Value', Dx);
        handles(17)=Dx;
        set(1,'Userdata', handles); 
        WPACKET('again');
       
elseif strcmp(Command, 'Run'),
	  runFlag=1;
      stopFlag = 0;
	  set(H_RUN, 'Enable', 'off');
      set(H_STOP, 'Enable', 'on');
      set(H_POP, 'Enable', 'off');
	  WPACKET('again');

elseif strcmp(Command, 'Stop'),
	  runFlag=0;
	  set(H_STOP, 'Enable', 'off');
      set(H_RUN, 'Enable', 'on');
      set(H_REST, 'Enable', 'on');
      set(H_POP, 'Enable', 'off');
      WPACKET('again');
      
elseif strcmp(Command, 'Rest'),
      runFlag=0;
      stopFlag =1;
      set(H_STOP, 'Enable', 'off');
      set(H_RUN, 'Enable', 'on');
      set(H_POP, 'Enable', 'on');
      set(H_REST, 'Enable', 'off');
      WPACKET('again');

elseif strcmp(Command, 'Res2Def'),
        Dx=0.9864;
        handles(17)=Dx;
        set(1,'Userdata', handles);    
        WPACKET('SetDx');     
        
elseif strcmp(Command, 'SetDx'),
        set(H_WED, 'String', num2str(Dx));  
        set(H_WSL, 'Value', Dx);
        Dx = handles(17);
        %Dx
        set(H_DxTxt, 'String', ['<\Delta x> = ', num2str(Dx),' ', 197]);
        set(H_WED, 'String', num2str(Dx));

        %[Psi, x, Ek, lambda] =  GSSPCKT(xo, Dx, Ko, relMass, xWidth, pNum);
         WPACKET('again');
    
elseif strcmp(Command, 'Choose_pot'),
      if runFlag,
         return;
      end   
      potIndex = get(H_POP, 'Value')-1;
      [V, xo, Dx, Comment]=PTNTLS(x, potIndex);
      handles(17) = Dx;
      set(1,'Userdata', handles);    
      WPACKET('SetDx');
      [Psi, x, Ek, lambda] =  GSSPCKT(xo, Dx, Ko, relMass, xWidth, pNum);
      runFlag=0;
      WPACKET('again');
      
elseif strcmp(Command, 'again'),
    	figure(1)
        
		if (runFlag==0)&(stopFlag==1), %#ok<ALIGN>
			[Psi, x, Ek, lambda] =  GSSPCKT(xo, Dx, Ko, relMass, xWidth, pNum);
			Eav = Ek+sum((abs(Psi).^2).*V);  % averaged state energy	 
			handles(20)=Eav;
            handles(21)=Ek;
            handles(22)=lambda;
			set(1, 'Userdata', handles);	
            set(H_LV, 'YData', V , 'Color', 'b', 'LineWidth', 2,'EraseMode', 'xor');
        end
				

        axes(H_AX);
          
		if abs(lambda)> 1e9,
			lambdaStr='\infty';
		else
			lambdaStr=[num2str(lambda),' ',197];
		end	
          
  
		% --------------- WAVE_PACKET ANIMATION --------------------------------	
		if Eav<50, %#ok<ALIGN>
			dt=0.01*4.13/50;
		else	
			dt = 0.01*4.13/real(Eav); % elementary time-interval in femtoseconds
			% it must be E*dt/h_bar << 2*pi, i.e. dt << 4.13/E femtosec 
        end
      
        set(H_TEav, 'Position', [-19, Eav+16], 'String', '<E>','Color', 'r');     
        set(H_LEav,'YData', [Eav, Eav]);
      
        set(H_EkTxt, 'String', ['<E_{kin}> = ',num2str(Ek),' eV']);
        set(H_EavTxt, 'String', ['<E> = ',num2str(Eav),' eV']); 
        set(H_LamTxt, 'String', ['<\lambda> = ',lambdaStr]);
        set(H_DxTxt, 'String', ['<\Delta x> = ', num2str(Dx),' ', 197]); 
      
        CycleNum=40; % number of times the equation is incrementally solved by the routine WVPCKT.MEX
		PAR=[dt, xWidth, relMass, CycleNum];
		j=0;
		firstPlot=1;	
        while runFlag==1 || firstPlot==1, %#ok<ALIGN>
            set(H_LR, 'YData', 400*imag(Psi));
            set(H_LI, 'YData', 400*real(Psi));             		
			set(H_TXT, 'String', num2str(j*dt*CycleNum));
            set(H_LP, 'YData',(0.05*400*400)*abs(Psi).^2);         
            drawnow;
			firstPlot=0;
            Psi=SCHROEQ(Psi, V, PAR);
            j=j+1;					
	   	end
end
