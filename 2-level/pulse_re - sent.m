%---Calculate number of cycles for a given pulse---
clear all
lambda = 800.;                 % in nm
lambda = lambda*1e-9;         % in m
bandw = 200;                % Pulse bandwidth in nm
bandw = bandw*1e-9;           % in m
N = 1024;                     % Number of points in the pulse
NN = 2^19;                    % Number of points in the Fourier result
c = 299792458;                % Speed of light in m/s
ph = 0;                       % Instanteneous pulse phase
wc = (2.*pi*c)/lambda;        % Center Frequency of Laser pulse
FWHM1 = 1e15*(lambda^2)/(c*bandw);  % Time FWHM in fs
offset = lambda - 2.*bandw;
start = 0.;
stop = 4.*bandw;
w = linspace(start,stop,N);
dw = w(2)-w(1);
w0 = (stop - start)/2.;
gaun1 = exp(-4*log(2.)*((w-w0)/bandw).^2);
gaun = gaun1.*(cos(wc*(w-w0)+ph));
ww = (w+offset)*1e9;                      % in nm (for plotting)
%---------------------------------------------------------
figure(1)
plot(ww, gaun, ww, gaun1);
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
%---------------------------------------------------------
gau = fftshift(fft(gaun1,NN));
t = (1:NN);   
[FWHM,xmin] = fwhm(t,gau);
peakpt = xmin + round(FWHM/2);
start = peakpt - 2*FWHM;
stop = peakpt + 2*FWHM;
gauf1 = abs(gau(start:stop));
Nf = length(gauf1);
tt = t(start:stop);
dt = FWHM1/FWHM;
tt = dt*(tt - tt(1));
tt0 = (tt(Nf)-tt(1))/2.;
gauf = gauf1.*(cos(wc*((tt-tt0)*1e-15)+ph));
%---------------------------------------------------------------
figure(2)
plot(tt, gauf, tt, gauf1), grid on;
xlabel('Time (fs)');
ylabel('Intensity (arbitrary units)');
title(sprintf('Pulse in Time Domain with FWHM = %d fs',FWHM1));
%---------------------------------------------------------------
nct = round((FWHM1*1e-15)*wc/(2.*pi));     % Number of cycles in time FWHM
display(FWHM1)
display(nct)