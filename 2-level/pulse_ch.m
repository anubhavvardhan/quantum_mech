%---Calculate number of cycles for a given pulse---
clear all
lambda = 800;                 % in nm
lambda = lambda*1e-9;         % in m
bandw = 21;                   % Pulse bandwidth in nm
bandw = bandw*1e-9;           % in m
N = 1000;                     % Number of points in the pulse (Length of 'w')
c = 299792458;                % Speed of light in m/s
ph = 0;                       % Instanteneous pulse phase
wc = (2.*pi*c)/lambda;        % Center Frequency of Laser pulse
FWHM1 = 1e15*(lambda^2)/(c*bandw);  % Time FWHM in fs
offset = lambda - 2.*bandw;
start = 0.;
stop = 4.*bandw;
w = linspace(start,stop,N);
dw = w(2)-w(1);                
Fs = dw/c;                    % Sampling Time period
w0 = (stop - start)/2.;
gaun = exp(-4*log(2.)*((w-w0)/bandw).^2).*exp(1i*(wc*(w-w0)+ph));
ww = (w+offset)*1e9;                      % in nm (for plotting)
oscgaun = length(crossing(real(gaun)));
display(oscgaun)
%---------------------------------------------------------
figure(3)
plot(ww, real(gaun), ww, imag(gaun), ww, abs(gaun));
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
%---------------------------------------------------------
% gau = fftshift(fft(gaun));
NN = 2^nextpow2(N); % Next power of 2 from length of 'w'
gau = ifft(gaun,NN)/N;
t = Fs/2*linspace(0,1,NN/2+1);
oscgau = length(crossing(real(gau)));
display(oscgau)
% t = (1:NN); 
% t = (1:N); 
[FWHM,xmin] = fwhm(t,gau);
% [FWHM,xmin] = findFWHM(t,gau);
display(FWHM)
peakpt = xmin + round(FWHM/2);
start = peakpt - 2*FWHM;
stop = peakpt + 2*FWHM;
gauf = gau(start:stop);
Nf = length(gauf);
tt = t(start:stop);
dt = FWHM1/FWHM;
tt = dt*(tt - tt(1));
%---------------------------------------------------------------
figure(4)
%plot(tt, real(gauf), tt, imag(gauf), tt, abs(gauf)), grid on;
plot(tt, real(gauf), tt, abs(gauf)), grid on;
xlabel('Time (fs)');
ylabel('Intensity (arbitrary units)');
title(sprintf('Pulse in Time Domain with FWHM = %d fs',FWHM1));
%---------------------------------------------------------------
nct = round((FWHM1*1e-15)*wc/(2.*pi));     % Number of cycles in time FWHM
display(FWHM1)
display(nct)
oscgauf = length(crossing(real(gauf)));
display(oscgauf)