clear all
bandwidth = 5;                    
bandwidth = bandwidth*1e-9;              % pulse bandwidth in m
lambda = 400.;                
lambda = lambda*1e-9;     	         % meter
N = 1024;                    		 % I saw on the net that it is usually taken as 1024, but would 512 be alright???
Fourierpoint = 2^19;          
hc = (6.62e-34)/(2.*pi);      
ph = 0;                           	 % Pulse phase
c = 299792458;                
fc = (2.*pi*c)/lambda;       		 % Center Frequency of Laser
d = 1;                      	         % dipole moment 
mrabi = 1;                    
FWHM1 = 1e15*(lambda^2)/(c*bandwidth);       % Time of Full Width at Half Maximum in fs
offset = lambda - 2.*bandwidth;
start = 0.;
stop = 4.*bandwidth;
w = linspace(start,stop,N);
dw = w(2)-w(1);
w0 = (stop - start)/2.;
% gaun1 = exp(-4*log(2.)*((w-w0)/bandwidth).^2);
% gaun = gaun1.*(cos(fc*(w-w0)+ph));
gaun = exp(-4*log(2.)*((w-w0)/bandwidth).^2).*exp(1i*(fc*(w-w0)+ph));
ww = (w+offset)*1e9;                      % in nm (for plotting)
%
%
figure(1)
% plot(ww, gaun, ww, gaun1);
plot(ww, real(gaun), ww, imag(gaun), ww, abs(gaun));
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
%
%
gau = fftshift(fft(gaun,Fourierpoint));
t = (1:Fourierpoint); 
[FWHM,xmin] = fwhm(t,gau);
peakpt = xmin + round(FWHM/2);
start = peakpt - 2*FWHM;
stop = peakpt + 2*FWHM;
gauf1 = abs(gau(start:stop));
% gauf = gau(start:stop);
Nf = length(gauf1);
tt = t(start:stop);
dt = 1e-15*FWHM1/FWHM;
% dt = FWHM1/FWHM;
tt = dt*(tt - tt(1));
tt0 = (tt(Nf)-tt(1))/2.;
gauf1=gauf1/max(gauf1);
gauf = gauf1.*(cos(fc*(tt-tt0)+ph));
%
%
figure(2)
plot(tt*1e15, gauf, tt*1e15, gauf1), grid on;
%plot(tt, gauf1), grid on;
%plot(tt, real(gauf), tt, imag(gauf), tt, abs(gauf)), grid on;
%plot(tt*1e15, real(gauf), tt*1e15, abs(gauf)), grid on;
% plot(tt, real(gauf), tt, abs(gauf)), grid on;
xlabel('Time (fs)');
ylabel('Intensity (arbitrary units)');
title(sprintf('Pulse in Time Domain with FWHM = %d fs',FWHM1));
nct = round((FWHM1*1e-15)*fc/(2.*pi));     % Number of cycles in time FWHM
display(FWHM1)
display(nct)
[re,rg] = evolution_RWA(Nf,gauf1,mrabi,fc,hc,dt);
figure(3)
plot(tt*1e15,abs(re),tt*1e15,abs(rg),tt*1e15, gauf1.^2)
% plot(tt*1e15,abs(re))
xlabel('Time (fs)');
ylabel('Population');
legend excited ground pulse^2
% legend excited
title(sprintf('Population Evolution in Time Domain with %d pi Pulse Area',mrabi));