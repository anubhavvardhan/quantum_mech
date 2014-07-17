	function POTWELLS(Command, radioButt)
%		 
%	The eigenvalue equation solved by the routine is:
%
%	 -K*[d^2u(x)]/dx^2 + V(x)*u(x) = u(x)*E 
%
%	with K = (h_bar)^2/(2*Me) =  3.8099827 eV*A^2 (A=Angstrom), where
%	h_bar is the Planck constant divided by 2*pi and Me is the electron mass.
%
%	Useful data:
%	h_bar= 0.65821 eV*fsec, the Planck's constant (fsec=femtosecond = 10^(-12) sec)
%	c = 2998 A/fsec, the velocity of light (A = Angstrom = 10^(-11) m)
%	m_e = 5.11e5 eV, the electron mass
%
%	The X-axis interval is discretised as a set of N=201 points. Initially,
%	the x-axisscale ranges from -10 A to 10 A (A=Angstrom)

global Psi Phi Pden E x N V DV Vo Vfig Vold runFlag stopFlag to


if nargin < 1,
   Command = 'Initialize';
	figure(1)
	clf
    set(1, 'Name', 'Potential wells -  R.Nobili - Padova University - ITALY');         
	set(1, 'Userdata', zeros(size(1:36)));
end

if ~strcmp(Command, 'Initialize'),
	% retrieving data saved at the end of the 'Initialize' call;    
	handle=get(1, 'Userdata');
   
   % retrieving data from the Userdata buffer of the window
   H_RAD=zeros(1,4);
   
   H_LSPEC=zeros(1:5);
   
   H_FIG= handle(1);
   H_AX = handle(2); 
   H_LV = handle(3);
   H_LP = handle(4);
   H_LE = handle(5);
   H_ED = handle(6);
   H_ETXT= handle(7);
   H_POP = handle(8);
   H_SLSCALE = handle(9);
   H_RAD(1)=handle(10);
	H_RAD(2)=handle(11);
 	H_RAD(3)=handle(12);
 	H_RAD(4)=handle(13);
   H_EDFREQ=handle(14); 
   H_STOP = handle(15);
   H_CONT = handle(16);
	H_REST = handle(17);
	H_LFR = handle(18);
 	H_LFI = handle(19);
   H_ETXT0 = handle(20);
   H_EDSCALE = handle(21);
   H_SLFREQ = handle(22);
	H_SPEC = handle(23); 
   
   dT = handle(24);
   % Vtop = handle(25);
   % Vbottom=handle(26);
   scale = handle(27);
   oldRadioButt=handle(28);
	N = handle(29);
	K = handle(30);
    M = handle(31);
	% oscFlag=handle(32);
   freq=handle(33);
   Yr = handle(34);
   level = handle(35);
	H_LSPEC(1) = handle(36);
   H_LSPEC(2) = handle(37);
   H_LSPEC(3) = handle(38);
   H_LSPEC(4) = handle(39);
   H_LSPEC(5) = handle(40);
	freqFlag = handle(41);
   H_FCBUTT = handle(42);
   H_FREQTXT = handle(43);
   
end
   
if strcmp(Command, 'Initialize'),
   
      % -------------- SYSTEM SETTING ------------------------------
      handle=zeros(1,46);
      
      K = 3.8099827;     
      M=80;
      level=1;
      scale = 0.1;
      x=scale*(-100:100);
      N=length(x);
      V=zeros(1,N);
      Vold=V;
      Vfig=V;
      runFlag=0;
      stopFlag=0; 
      freqFlag=1; % frequency change enabled 
      oscFlag=0;
      to=0;
      dT=0.001;      
      freq=10;
      Phi = zeros(1,N);
      Psi = zeros(N, N);
      Pden =zeros(N, N);
      
      Vbottom=-2;
      Vtop=4;
      Yr=0;
      
      H_AX = axes('Units', 'normalized',...
         'Position', [0.075, 0.15, 0.725, 0.78], ...
			'XLim',[x(1),x(N)], ...
         'YLim',[Vbottom, Vtop],...
         'Box', 'on');
      
      
      title('Energy eigenfunctions and eigenvalues of an electron');
   	  xlabel(['Electron position (Angstr',246,'m)']);
	  ylabel('Energy level (eV)');   
      

      H_LV = line('XData', x, 'Ydata', V, 'Color',...
   				'b', 'LineWidth', 2, 'EraseMode', 'xor');
         		            
      H_LP = line('XData', x, 'Color','k', 'EraseMode', 'xor');
      
      H_LFR = line('XData', x, 'Color',[1,0.75, 1],'EraseMode', 'xor');		

	  H_LFI = line('XData', x, 'Color','c','EraseMode', 'xor');
		      
      H_LE = line('XData', [x(1), x(N)], 'Color','r', 'EraseMode', 'xor');
  
  		line('XData', [100*x(1), 100*x(N)],... % baseline of probability density
              'YData', [0, 0],...
              'LineStyle', ':',...
     			  'Color',[0.8, 0.8, 0.8], 'EraseMode', 'xor');
             
      if Yr ~=0,      %#ok<ALIGN>
         line('XData', [100*x(1), 100*x(N)],... % baseline of complex and real parts
              'YData', [Yr, Yr],...
              'LineStyle', ':',...
     			  'Color', [0.8, 0.8, 0.8], 'EraseMode', 'xor');
		end
  	  
  % ---------------- UICONTROLS -------------
  
    
  	  H_FIG=1;
      uicontrol(H_FIG, 'style', 'text',...
         'string', 'Energy level No:', ...
         'Units', 'normalized',...
         'Position', [0.7, 0.015, 0.175, 0.04],...
         'Background', [0.8, 0.8, 0.8]);
      
      text(0.75, 0.15, '- Probability density',...
         'Unit', 'Normalized',...
         'Color','k', 'FontName', 'Arial',...
         'FontSize', 8);
      
      
      text(0.75, 0.1, '- Real part',...
         'Unit', 'Normalized',...
      	'Color','m', 'FontName', 'Arial',...
         'FontSize', 8);

		text(0.75, 0.05, '- Imaginary part',...
   		'Unit', 'Normalized',...
      	'Color',[0, 0.8, 0.8], 'FontName', 'Arial',...
         'FontSize', 8);

	
      H_ED =  uicontrol(H_FIG, 'style', 'edit',...
         'Units', 'normalized',...
         'Position', [0.875, 0.02, 0.05, 0.04],...
         'String', '1',...
      	'Background', [0.8, 0.8, 0.8]);
      
      %H_BUP =  
         uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'up',...
         'Callback', 'POTWELLS(''UpButt'');',...
         'Units', 'normalized',...
         'Position', [0.93, 0.04, 0.06, 0.04]);
      
      % H_BDOWN =  
          uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'down',...
         'Callback', 'POTWELLS(''DownButt'');',...
         'Units', 'normalized',...
         'Position', [0.93, 0.0, 0.06, 0.04]);

 		H_ETXT0= uicontrol(H_FIG, 'style', 'text',...
         'string', 'Energy (eV) = ', ...
         'HorizontalAlignment', 'right',...
         'Units', 'normalized',...
         'Position', [0.12, 0.155, 0.15, 0.03],...
         'Background', 'w'); %[0.8, 0.8, 0.8]);
      
		H_ETXT = uicontrol(H_FIG, 'style', 'text',...
         'string', '0.0', ...
         'HorizontalAlignment', 'left',...
         'Units', 'normalized',...
         'Position', [0.28, 0.155, 0.15, 0.03],...
         'Background', 'w'); %[0.8, 0.8, 0.8]);
      
		Yo=0.105;	      

        uicontrol(H_FIG, 'style', 'text',...
         'string','Potential type', ...
         'Units', 'normalized',...
         'Position', [0.8, 0.865+Yo, 0.18, 0.03],...
         'Background', [0.8, 0.8, 0.8]);

      
		H_POP = uicontrol(H_FIG, 'style', 'pop',...
         'string',['Flat well|Step well|Double well|Triple well|'...
         'Flat well tilted|Step well tilted|Double well tilted|'...   
         'Staircase|Staircase tilted|Multibarrier|'...
         'Parabolic|Quartic|Quartic tilted|Parabolic in expansion'],...
         'Callback', 'POTWELLS(''Choose_pot'')',...
         'Value', 1,...
      	 'Units', 'normalized',...
         'Position', [0.8, 0.84+Yo, 0.2, 0.02],...
         'Background', 'w');
      
	      
        H_SLSCALE = uicontrol(H_FIG, 'style', 'slider',...
            'Min', 0.01, 'Max', 0.8, 'Value', scale, ...
            'Callback','POTWELLS(''ScaleSlider'')', ...
            'Units', 'normalized',...
            'Position', [0.41, 0.01, 0.14, 0.04]);
       
		H_EDSCALE = uicontrol(H_FIG, 'style', 'edit',...
            'Callback','POTWELLS(''EdScale'')', ...
 			'Units', 'normalized',...
            'Position', [0.35, 0.01, 0.06, 0.04],...
            'String', num2str((N-1)*scale),... 
            'Background', 'w');
      
   	uicontrol(H_FIG, 'style', 'text',...
         'string','X-axis width (Angstrom):', ...
         'Units', 'normalized',...
         'Position', [0.14, 0.01, 0.21, 0.04],...
         'Background', [0.8, 0.8, 0.8]);

   
      % ---------- DYNAMIC SECTION ---------------------
	   uicontrol(H_FIG, 'style', 'frame',...
         'Units', 'normalized',...
         'Position', [0.805, 0.31+Yo, 0.195, 0.495],...
         'Background', [0.8, 0.8, 0.8]);

      
          
		uicontrol(H_FIG, 'style', 'text',...
         'string','DYNAMIC SECTION', ...
         'Units', 'normalized',...
         'Position', [0.81, 0.773+Yo, 0.18, 0.03],...
         'Background', [0.8, 0.8, 0.8]);
      
  		H_RAD=zeros(1,4);

 		H_RAD(1) = uicontrol(H_FIG, 'style', 'radio',...
         'string','Fast change', ... %         'Value', 1, ...
      	'Callback','POTWELLS(''Radio'', 1)', ...
      	'Units', 'normalized',...
      	'Position', [0.81, 0.73+Yo, 0.18, 0.03],...
         'Background', [0.8, 0.8, 0.8]);
      
        H_RAD(2) = uicontrol(H_FIG, 'style', 'radio',...
            'string','Moderate change', ...
            'Callback','POTWELLS(''Radio'', 2)', ...
            'Units', 'normalized',...
            'Position', [0.81, 0.68+Yo, 0.185, 0.03],...
            'Background', [0.8, 0.8, 0.8]);
      
      H_RAD(3) = uicontrol(H_FIG, 'style', 'radio',...
         'string','Slow change', ...
         'Callback','POTWELLS(''Radio'', 3)', ...
	     'Units', 'normalized',...
         'Position', [0.81, 0.63+Yo, 0.185, 0.03],...
         'Background', [0.8, 0.8, 0.8]);
      
 	   H_RAD(4) = uicontrol(H_FIG, 'style', 'radio',...
         'string','10 oscillations', ...
         'Callback','POTWELLS(''Radio'', 4)', ...
         'Units', 'normalized',...
         'Position', [0.81, 0.58+Yo, 0.18, 0.03],...
         'Background', [0.8, 0.8, 0.8]);
      
       H_SLFREQ = uicontrol(H_FIG, 'style', 'slider',...
         'Min', 0.01, 'Max', 10, 'Value', 0.05, ...
         'Callback','POTWELLS(''FreqSlider'')', ...
 		 'Units', 'normalized',...
         'Position', [0.87, 0.5+Yo, 0.125, 0.04]);
      
      H_EDFREQ = uicontrol(H_FIG, 'style', 'edit',...
         'Callback','POTWELLS(''EdFreq'')', ...
 		 'Units', 'normalized',...
         'Position', [0.81, 0.5+Yo, 0.06, 0.04],...
         'String', ' ',... 
         'Background', 'w');
      
      H_FREQTXT=uicontrol(H_FIG, 'style', 'text',...
         'Units', 'normalized',...
         'Position', [0.81, 0.45+Yo, 0.185, 0.04],...
         'String', ' ',... 
         'Background', [0.8, 0.8, 0.8]);
      
      H_FCBUTT = uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'Enable freq. change',...
         'Callback', 'POTWELLS(''FreqChange'');',...
         'Units', 'normalized',...
         'Position', [0.81, 0.41+Yo, 0.185, 0.04]);

      uicontrol(H_FIG, 'style', 'text',...
         'string','Frequency (1/fsec)', ...
         'Units', 'normalized',...
         'Position', [0.81, 0.545+Yo, 0.185, 0.03],...
         'HorizontalAlignment', 'left',...
      	 'Background', [0.8, 0.8, 0.8]); 

      
      H_STOP =  uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'Stop',...
         'Enable', 'off',...
      	 'Callback', 'POTWELLS(''StopButt'');',...
         'Units', 'normalized',...
         'Position', [0.81, 0.36+Yo, 0.09, 0.04]);
      
      H_CONT =  uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'Continue',...
         'Enable', 'off',...
         'Callback', 'POTWELLS(''ContButt'');',...
         'Units', 'normalized',...
         'Position', [0.905, 0.36+Yo, 0.09, 0.04]);
      
      H_REST =  uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'Restore',...
         'Enable', 'off',...
         'Callback', 'POTWELLS(''RestoreButt'');',...
         'Units', 'normalized',...
         'Position', [0.81, 0.315+Yo, 0.09, 0.04]);
      
      H_SPEC =  uicontrol(H_FIG, 'Style', 'pushbutton',...
         'String', 'Spectrum',...
         'Enable', 'off',...
         'Callback', 'POTWELLS(''SpecButt'');',...
         'Units', 'normalized',...
         'Position', [0.905, 0.315+Yo, 0.09, 0.04]);     
      
                
      % H_ASPEC = 
       axes('Units', 'normalized',...
         'Position', [0.805, 0.15+Yo, 0.195, 0.162], ...
			'XLim',[0.5,5.5], ...
         'YLim',[0,1],...
         'Box', 'on',...
         'YTickMode', 'manual',...
         'FontName', 'Arial',...
         'FontSize', 7,...
         'XTickLabel', 'E1|E2|E3|E4|E5');
      
  		uicontrol(H_FIG, 'style', 'text',...
         'string','Energy spectrum of the system in the final potential configuration.', ...
         'Units', 'normalized',...
         'Position', [0.82, 0.037+Yo, 0.19, 0.084],...
         'HorizontalAlignment', 'left',...
         'Background', [0.8, 0.8, 0.8],...
         'FontName', 'Arial',...
         'FontSize', 7);

     	for n=1:5,       %#ok<ALIGN>
            H_LSPEC(n)= line('XData', [n,n], ...
            'YData', [0,0], 'Color', 'r');
		end
			     
      
      oldRadioButt=0;
      % time=0;
        handle(1) = H_FIG;
		handle(2) = H_AX; 
        handle(3) = H_LV;
        handle(4) = H_LP;
        handle(5) = H_LE;
        handle(6) = H_ED;
        handle(7) = H_ETXT;
        handle(8)= H_POP;
        handle(9) = H_SLSCALE;
        handle(10) = H_RAD(1);
 		handle(11) = H_RAD(2);
 		handle(12) = H_RAD(3);
        handle(13) = H_RAD(4);
        handle(14) = H_EDFREQ; 
 		handle(15) = H_STOP;
 		handle(16) = H_CONT;
 		handle(17) = H_REST;
		handle(18) = H_LFR;
        handle(19) = H_LFI;
		handle(20) = H_ETXT0;
        handle(21) = H_EDSCALE;
        handle(22) = H_SLFREQ;
		handle(23) = H_SPEC;
        handle(24) = dT;
        handle(25) = Vtop;
        handle(26) = Vbottom;
        handle(27) = scale;
        handle(28) = oldRadioButt;
        handle(29) = N;		
        handle(30)= K;
        handle(31) = M;
 		handle(32) = oscFlag;
		handle(33) = freq;
        handle(34)=Yr;
        handle(35) = level;
        handle(36) = H_LSPEC(1);
        handle(37) = H_LSPEC(2);
        handle(38) = H_LSPEC(3);
        handle(39) = H_LSPEC(4);
        handle(40) = H_LSPEC(5);
        handle(41) = freqFlag;
		handle(42) = H_FCBUTT;
        handle(43) = H_FREQTXT;
      

        set(H_FIG, 'Userdata', handle); % save data in Userdata buffer
      
        POTWELLS('Choose_pot');
      
   elseif strcmp(Command, 'EdScale'),
      if runFlag,
         set(H_EDSCALE, 'String', num2str(scale));
         return
      end   
      scaleString = get(H_EDSCALE, 'String');
      scale = str2double(scaleString);
      maxVal=get(H_SLSCALE, 'Max');
      minVal=get(H_SLSCALE, 'Min');
      if scale > maxVal, %#ok<ALIGN>
         scale = maxVal;
         set(H_EDSCALE, 'String', maxVal);
      elseif scale < minVal,
			scale = minVal;
         set(H_EDSCALE, 'String', minVal);
		end

      set(H_SLSCALE, 'Value', scale);
      handle(27)=scale;
      set(H_FIG, 'Userdata', handle);
      POTWELLS('ScaleSlider');
      
   elseif strcmp(Command, 'ScaleSlider'),
      if runFlag,
         set(H_SLSCALE, 'Value', scale);
         return
      end   
      scale = get(H_SLSCALE, 'Value');
      set(H_EDSCALE, 'String', num2str((N-1)*scale));

      handle(27)=scale;
      
      set(1, 'Userdata', handle); 
      
      x = scale*(-100:100);

      set(H_AX, 'XLim', [x(1), x(N)]);
      set(H_LV, 'Xdata', x);
      set(H_LP, 'Xdata', x);
      set(H_LFR,'Xdata', x);
      set(H_LFI,'Xdata', x);
      set(H_LE, 'Xdata', [x(1), x(N)]);

      POTWELLS('Choose_pot');
           
    
   elseif strcmp(Command, 'EdFreq'),
      if runFlag,
         set(H_EDFREQ, 'String', num2str(freq));
         return
      end         
      freqString = get(H_EDFREQ, 'String');
	  freq = str2double(freqString);
      maxVal=get(H_SLFREQ, 'Max');
      minVal=get(H_SLFREQ, 'Min');
      if freq > maxVal, %#ok<ALIGN>
         freq = maxVal;
         set(H_EDFREQ, 'String', maxVal);
      elseif freq < minVal,
			freq = minVal;
         set(H_EDFREQ, 'String', minVal);
		end

      set(H_SLFREQ, 'Value', freq);
      handle(33)=freq;
      set(H_FIG, 'Userdata', handle);
      
 	elseif strcmp(Command, 'FreqSlider'),
      if runFlag,
         set(H_SLFREQ, 'Value', freq);
         return
      end   
      freq = get(H_SLFREQ, 'Value');
      set(H_EDFREQ, 'String', num2str(freq));

      handle(33)=freq;
      
      set(1, 'Userdata', handle); 
      
	elseif strcmp(Command, 'FreqChange'),
      if runFlag, %#ok<ALIGN>
            return
      end
      set(H_FREQTXT, 'String','');   
      if freqFlag==1,
         set(H_FCBUTT, 'String', 'Disable freq. change'); 
         freqFlag = 0; 
      else
         set(H_FCBUTT, 'String', 'Enable freq. change');
         freqFlag = 1;
      end   
      handle(41)=freqFlag;
      set(1, 'Userdata', handle); 

   elseif strcmp(Command, 'UpButt'),
      if runFlag && stopFlag,
         return
      end   
      level=level+1;
      if level > M,
         level = M;
      end   
      handle(35)=level;
      set(1, 'Userdata', handle);
      POTWELLS('Draw');
      
   elseif strcmp(Command, 'DownButt'),   
      if runFlag && stopFlag,
         return
      end   

      level=level-1;
      if level < 1,
         level = 1;
      end  
      
      handle(35)=level;
      set(1, 'Userdata', handle);
      POTWELLS('Draw');
      
   elseif strcmp(Command, 'Choose_pot'),
      if runFlag,
         return
      end   
      
      choice = get(H_POP, 'Value');
      % string = get(H_POP, 'String');
      
      Vo=zeros(1,N);
      V=Vo; 
      
         % 'Flat well|Step well|Double well|Triple well|'...
         % 'Flat well tilted|Step well tilted|Double well tilted|'...   
         % 'Staircase|Staircase tilted|Multibarrier|'...
         % 'Parabolic|Quartic|Quartic tilted|'...
       
      if choice==1, % FLAT WELL   
		       
			DV=(0.01/scale)*x;
         
      elseif choice==2, % STEP WELL    
         Na = floor(N/4);
    		V(Na:N)=ones(size(Vo(Na:N)));
			DV=(0.01/scale)*x;
            
      elseif choice==3, % DOUBLE WELL
         
         Na = floor(N/3);
         Nb = N-Na;  
         V(Na:Nb)= 2*ones(size(Vo(Na:Nb)));           
         V(Nb:N)=ones(size(Vo(Nb:N)));
      	 DV=(0.01/scale)*x;
   
      elseif choice==4, % TRIPLE WELL
         
         Na = floor(N/4);
         Nb = floor(N/2);
         Nc = N-Na;  
         V(Na:Nb)= 2*ones(size(Vo(Na:Nb)));           
         V(Nb:Nc)=-1*ones(size(Vo(Nb:Nc)));
         V(Nc:N)= ones(size(Vo(Nc:N)));
			DV=(0.01/scale)*x;
         
      elseif choice==5, % FLAT WELL TILTED
         V = V - (0.01/scale)*x;
			DV=(0.02/scale)*x;
             
      elseif choice==6, % STEP WELL TILTED
         
         Na = floor(N/4);         
		 V(Na:N)=ones(size(Vo(Na:N)));
       	 V = V - (0.01/scale)*x;
         DV=(0.02/scale)*x;
 
       elseif choice==7, % DOUBLE WELL TILTED
         
         Na = floor(N/3);
         Nb = N-Na;  
         V(Na:Nb)= 2*ones(size(Vo(Na:Nb)));           
         V(Nb:N)=ones(size(Vo(Nb:N)));
         V = V - (0.01/scale)*x;
 		 DV=(0.02/scale)*x;
        
       elseif choice==8, % STAIRCASE 
         Na = floor(N/6);
		 Nb = floor(2*N/6);
		 Nc = floor(3*N/6);
         Nd = floor(4*N/6);
         V(Na:Nb)= 0.5*ones(size(Vo(Na:Nb)));
         V(Nb:Nc)= 1.5*ones(size(Vo(Nb:Nc)));
		 V(Nc:Nd)= 2*ones(size(Vo(Nc:Nd)));
		 V(Nd:N)= 2.5*ones(size(Vo(Nd:N)));
      	 DV=(0.01/scale)*x;
	              
      elseif choice==9, % STAIRCASE TILTED
         
         Na = floor(N/6);
		 Nb = floor(2*N/6);
		 Nc = floor(3*N/6);
         Nd = floor(4*N/6);
         V(Na:Nb)= 0.5*ones(size(Vo(Na:Nb)));
         V(Nb:Nc)= 1.5*ones(size(Vo(Nb:Nc)));
		 V(Nc:Nd)= 2*ones(size(Vo(Nc:Nd)));
		 V(Nd:N)= 2.5*ones(size(Vo(Nd:N)));
		 V = V - (0.01/scale)*x;
      	 DV=(0.02/scale)*x;
	   
      elseif choice==10, % MULTIBARRIER
            V=Vo;
            
         	W1 =floor(N/21);
			W2=floor(N/33);
         	Nz=zeros(1,26);
         
         	Nz(1)=W1-1;
         	Nz(2)=Nz(1)+W2;
         	V(Nz(1):Nz(2))= ones(size(Vo(Nz(1):Nz(2))));
         
         	for r=3:2:25,
            	Nz(r)=Nz(r-1)+W1;
            	Nz(r+1)=Nz(r)+W2;
				V(Nz(r):Nz(r+1))= ones(size(Vo(Nz(r):Nz(r+1))));
         	end         
				DV=(0.01/scale)*x;
  	
        elseif choice==11, % PARABOLIC

			V = x.^2;
      	    V=(4/V(1))*V;
			DV=(0.01/scale)*x;
         
        elseif choice==12, % QUARTIC
         
         const = x(1)^4;
         fact= 0.8*x(1)*x(1);
         V = x.^4 - fact*x.^2 + 0.1*const;
         V=(4/V(1))*V;
			DV=(0.01/scale)*x;
         
            
       elseif choice==13, % QUARTIC TILTED   
         
         const = x(1)^4;
         fact= 0.8*x(1)*x(1);
         V = x.^4 - fact*x.^2 + 0.1*const;
         V=(4/V(1))*V;
         V = V - (0.01/scale)*x;
		 V=(4/V(1))*V;
         DV=(0.02/scale)*x;
            
       elseif choice==14, % QUARTIC TILTED   
         
           	V = x.^2;
      	    V=(32/V(1))*V;
			DV=(0.01/scale)*x;
       
       end   
      
      V(1)=1e6;
      V(N)=1e6;
      
     
      set(H_LV, 'Ydata', V);
      
      POTWELLS('Solve_equation');
      POTWELLS('Draw');
      
   elseif strcmp(Command, 'Solve_equation'),
      
      if runFlag,
         return
      end   

		Laplacian = diag(ones(1,N-1),1)+diag(ones(1,N-1),-1)-diag(2*ones(1,N));

      H = -K/(scale*scale)*Laplacian+diag(V);

		[Psi,D]=eig(H);
		[E, Ind]=sort(diag(D));

		Psi=Psi(:,Ind)+i*eps*zeros(size(Psi));
          
      % D=D(:,Ind);
      
      M=40;

		for n=1:M, %#ok<ALIGN>
   		Psi(:,n)=Psi(:,n)/norm(Psi(:,n), 'fro');
   		a=sum(Psi(1:floor(N/2), n));	
   		if a < 0,
      		Psi(:, n)=-Psi(:,n);
   		end   
      	Pden(:, n) = abs(Psi(:,n)).^2;   
      end      
      
		            
     elseif strcmp(Command, 'Draw'),
	  if runFlag, %#ok<ALIGN>
         return
      end   
            
      set(H_LP, 'Ydata', 60*Pden(:,level),'EraseMode', 'xor');
      
      set(H_LFR,'Ydata', Yr+4*real(Psi(:, level)));
            
      set(H_LFI,'Ydata', Yr+4*imag(Psi(:, level)));
      
      
      set(H_LE, 'Ydata', [E(level), E(level)], 'Color','r', 'EraseMode', 'xor');      
      
      set(H_ED, 'String', num2str(level));
      
      set(H_ETXT, 'String', num2str(E(level)));
      
   elseif strcmp(Command, 'Radio'),
      if runFlag,
      	POTWELLS('RestoreButt');
      end
      h_bar= 0.65821; % eV*fsec
      
      set(H_FREQTXT, 'String','');
      if freqFlag,        %#ok<ALIGN>
      	if radioButt==1, %#ok<ALIGN>
            freq = 1e-2/dT;  % help schrstp.m
         elseif radioButt==2,
         	freq=1e-4/dT;
      	elseif radioButt==3,
         	freq=1e-5/dT;
      	elseif radioButt==4,
         	freq=(E(level+1)-E(level))/(2*pi*h_bar); % help schrosc.m
         	set(H_FREQTXT, 'String',['Freq.=[E(', num2str(level+1),')-E(',num2str(level),')]/h']) 
   
         end
		end

      set(H_RAD, 'Value', 0);
 		set(H_RAD(radioButt), 'Value', 1);
	   handle(28)=radioButt;
      handle(33)=freq;
      
      set(H_FIG, 'Userdata', handle);
      
      set(H_SLFREQ, 'Value', freq);
      set(H_EDFREQ, 'String', num2str(freq));

      
      runFlag=1;
      stopFlag=0;
      % time=0;
    
      set(H_LE, 'Ydata', [10, 10]);         
		set(H_ETXT0, 'String', 'Time (fsec) = ');
	
		set(H_STOP, 'Enable', 'on');      
      set(H_FCBUTT, 'Enable', 'off');      
      
      POTWELLS('Solve_equation');
	  POTWELLS('InitDynamics', radioButt);
      POTWELLS('Evolve', radioButt);
     
   elseif strcmp(Command, 'StopButt'),
     
      set(H_RAD, 'Value', 0);
      runFlag=0;
      stopFlag=1;
      
      set(H_CONT, 'Enable', 'on');
      set(H_REST, 'Enable', 'on');
      set(H_SPEC, 'Enable', 'on');
      set(H_STOP, 'Enable', 'off');
      
      for n=1:5,
         set(H_LSPEC(n), 'YData', [-1,-1]);
      end   

   elseif strcmp(Command, 'ContButt'),
      
      for n=1:5,
         set(H_LSPEC(n), 'YData', [-1,-1]);
      end  

	  set(H_SPEC, 'Enable', 'off');      
      set(H_STOP, 'Enable', 'on');
      set(H_REST, 'Enable', 'off');
  	  set(H_SPEC, 'Enable', 'off');

      set(H_RAD(oldRadioButt), 'Value', 1);
      runFlag=1;	 
      
      POTWELLS('Evolve', oldRadioButt);
      
   elseif strcmp(Command, 'RestoreButt'),
      
      runFlag=0;	
      stopFlag=0;
      to=0;
      set(H_STOP, 'Enable', 'off');
      set(H_CONT, 'Enable', 'off');
	  set(H_REST, 'Enable', 'off');
      set(H_SPEC, 'Enable', 'off');
	  set(H_FCBUTT, 'Enable', 'on');      
      
      for n=1:5,
         set(H_LSPEC(n), 'YData', [-1,-1]);
      end  
      
      set(H_ETXT0, 'String', 'Energy (eV) = ');
      set(H_ETXT, 'String','0.00');

      set(H_RAD(oldRadioButt), 'Value', 0);
      
      set(H_LV, 'Ydata', V);
      
      POTWELLS('Solve_equation');
      POTWELLS('Draw');			
      
   elseif strcmp(Command, 'SpecButt'),
       Vold=V;   
      % V=Vfig;
      POTWELLS('Solve_equation');
      
      W=zeros(1:5);
      
      
      for n=1:5,
         W(n)=abs(Psi(:,n)'*conj(Phi))^2;      
      	 set(H_LSPEC(n),'YData',[0, W(n)]);      
      end
              
		V=Vold;
		POTWELLS('Solve_equation');
            
 	elseif strcmp(Command, 'InitDynamics'),
       
       
      Phi = Psi(:, level);
		         
   elseif strcmp(Command, 'Evolve'),
      
        
      if radioButt==4,
         cycleNum=100;
      else
         cycleNum=100;
      end
      
      omega =2*pi*freq; 
      tLim=10/freq;    
      
      PAR=[dT, omega, cycleNum, x(N)-x(1), 1];
      Zr=zeros(size(DV)); 
      
      while runFlag, %#ok<ALIGN>
      
      if radioButt==4, %#ok<ALIGN>
         if to<tLim,
            [Phi, to, Vfig] = Schrosc(Phi, to, V, DV, PAR);  
         else
            [Phi, to, Vfig] = schrosc(Phi, to, V, Zr, PAR);
         end   
      else
         [Phi, to, Vfig] =schrstp(Phi, to, V, DV, PAR);
		end
            
       set(H_LV, 'Ydata', Vfig);	     
       set(H_LFI,'Ydata', Yr + 4*imag(Phi));
       set(H_LFR,'Ydata', Yr + 4*real(Phi));   
       set(H_LP,'Ydata', 60*abs(Phi).^2);
 		 set(H_ETXT, 'String',num2Str(to));
	
       drawnow;
    end
end   
