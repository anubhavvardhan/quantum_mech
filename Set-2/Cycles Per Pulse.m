clear all
lambda = 800;                 
lambda = lambda*1e-9;         
bandw = 21;                   
bandw = bandw*1e-9;          
N = 1000;                     
c = 299792458;                
phase = 0;                       
wc = (2.*pi*c)/lambda;        
FWHM1 = 1e15*(lambda^2)/(c*bandw);  
offset = lambda - 2.*bandw;
start = 0.;
stop = 4.*bandw;
w = linspace(start,stop,N);
dw = w(2)-w(1);                
Fs = dw/c;                    
w0 = (stop - start)/2.;
gaun = exp(-4*log(2.)*((w-w0)/bandw).^2).*exp(1i*(wc*(w-w0)+phase));
ww = (w+offset)*1e9;                      
oscgaun = length(crossing(real(gaun)));
display(oscgaun)
%
figure(3)
plot(ww, real(gaun), ww, imag(gaun), ww, abs(gaun));
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
%
NN = 2^nextpow2(N); 
gau = ifft(gaun,NN)/N;
t = Fs/2*linspace(0,1,NN/2+1);
oscgau = length(crossing(real(gau)));
display(oscgau)
[FWHM,xmin] = fwhm(t,gau);
% 
display(FWHM)
peakpt = xmin + round(FWHM/2);
start = peakpt - 2*FWHM;
stop = peakpt + 2*FWHM;
gauf = gau(start:stop);
Nf = length(gauf);
tt = t(start:stop);
dt = FWHM1/FWHM;
tt = dt*(tt - tt(1));
%
figure(4)
plot(tt, real(gauf), tt, abs(gauf)), grid on;
xlabel('Time (fs)');
ylabel('Intensity (arbitrary units)');
title(sprintf('Pulse in Time Domain with FWHM = %d fs',FWHM1));
%
cycles = round((FWHM1*1e-15)*wc/(2.*pi));     
display(FWHM1)
display(cycles)
oscgauf = length(crossing(real(gauf)));
display(oscgauf)