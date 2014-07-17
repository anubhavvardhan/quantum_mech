fs = 10000;                          % Sample frequency (Hz)
% t = 0:1/fs:10-1/fs;                % 10 sec sample
% x = (1.3)*sin(2*pi*15*t) ...       % 15 Hz component
%   + (1.7)*sin(2*pi*40*(t-2)) ...   % 40 Hz component
%   + (2.5)*randn(size(t));          % Gaussian noise;

FWHM=0.5;
start=0;
stop=4*FWHM;
N=(stop-start)*fs;
t=linspace(start,stop,N);
t0=(stop-start)/2;
x=exp(-4*log(2)/(FWHM^2)*(t-t0).^2).*exp(100.*1i*(t-t0));
figure(1)
plot(t,real(x))
xlabel('Time (s)')
ylabel('Intensity')
title('{\bf Time-Instensity Plot}')
hold on
m = length(x);          % Window length
n = pow2(nextpow2(m));  % Transform length
y = fftshift(fft(x,n)); % DFT
f = (0:n-1)*(fs/n);     % Frequency range
power = y.*conj(y)/n;   % Power of the DFT
figure(2)
plot(f,power,f,real(y))
xlabel('Frequency (Hz)')
ylabel('Power')
title('{\bf Periodogram}')
hold on