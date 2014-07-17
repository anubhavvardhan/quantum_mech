clear all
c = 299792458;                
Fs = 1000;                    
T = 1/Fs;                     
L = 1000;                     
t = (0:L-1)*T;                
t = t - t(L/2);
y = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t); 
figure(1)
plot(Fs*t,y)
title('Signal Corrupted with Zero-Mean Random Noise')
xlabel('time (milliseconds)')
NFT = 2^nextpow2(L); 
Y = fft(y,NFT)/L;
f = Fs/2*linspace(0,1,NFT/2+1);

% Single-sided amplitude spectrum.
figure(2)
plot(f,2*abs(Y(1:NFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')