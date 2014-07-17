% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.

x = 1:NN;
y = gaun;
y = y / max(y);
N = length(y);
lev = 0.5;
if y(1) < lev                  % index of center of pulse
    [garbage,centerindex]=max(y);
    Pol = +1;
    disp('Pulse Polarity = Positive')
else
    [garbage,centerindex]=min(y);
    Pol = -1;
    disp('Pulse Polarity = Negative')
end
i = 2;
while sign(y(i)-lev) == sign(y(i-1)-lev)
    i = i+1;
end                                   %first crossing
crossing = (lev-y(i-1)) / (y(i)-y(i-1));
lead = x(i-1) + crossing*(x(i)-x(i-1));
i = centerindex+1;                    %next crossing at center
while ((sign(y(i)-lev) == sign(y(i-1)-lev)) && (i <= N-1))
    i = i+1;
end
if i ~= N
    Ptype = 1;  
    disp('Pulse is Impulse or Rectangular with 2 edges')
    crossing = (lev-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + crossing*(x(i)-x(i-1));
    width = ttrail - lead;
else
    Ptype = 2; 
    disp('Step-Like Pulse, no second edge')
    ttrail = NaN;
    width = NaN;
end
disp(width)