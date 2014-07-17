clear all
fs = 1e+10;                          % Sample frequency (Hz)
% t = 0:1/fs:10-1/fs;                % 10 sec sample
% x = (1.3)*sin(2*pi*15*t) ...       % 15 Hz component
%   + (1.7)*sin(2*pi*40*(t-2)) ...   % 40 Hz component
%   + (2.5)*randn(size(t));          % Gaussian noise;
lambda = 800;                 % in nm
lambda = lambda*1e-9;         % in m
bandw=21;
bandw = bandw*1e-9;           % in m
c = 299792458;                % Speed of light in m/s
ph = 0;                       % Instanteneous pulse phase
wc = (2.*pi*c)/lambda;        % Center Frequency of Laser pulse
FWHM1 = 1e15*(lambda^2)/(c*bandw);  % Time FWHM in fs
offset = lambda - 2.*bandw;
start=0;
stop=4*bandw;
N=(stop-start)*fs;
w=linspace(start,stop,N);
w0=(stop-start)/2;
gau=exp(-4*log(2.)*((w-w0)/bandw).^2).*exp(1i*(wc*(w-w0)+ph));
ww = (w+offset)*1e9;                      % in nm (for plotting)
figure(1)
plot(ww,real(gau),ww,abs(gau))
xlabel('Wavelength (nm)')
ylabel('Power')
title('{\bf Periodogram}')
hold on
m = length(gau);          % Window length
n = pow2(nextpow2(m));  % Transform length
y = n*ifftshift(ifft(gau,n)); % DFT
t = (0:n-1)*(fs/n);     % Frequency range
power = y.*conj(y)/n;   % Power of the DFT
figure(2)
plot(t,power,t,real(y))
xlabel('Time (s)')
ylabel('Intensity')
title('{\bf Time-Instensity Plot}')
hold on