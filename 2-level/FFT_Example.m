clear all
c = 299792458;                % Speed of light in m/s
Fs = 1000;                    % Sampling frequency
T = 1/Fs;                     % Sample time
L = 1000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
t = t - t(L/2);
y = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
figure(1)
plot(Fs*t,y)
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure(2)
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')