clear all
fs = 1e+10;                          % frequency
% t = 0:1/fs:10-1/fs;                
% x = (1.3)*sin(2*pi*15*t) ...       
%   + (1.7)*sin(2*pi*40*(t-2)) ...   
%   + (2.5)*randn(size(t));          
lambda = 800;                 
lambda = lambda*1e-9;         
bandw=21;
bandw = bandw*1e-9;           
c = 299792458;                
ph = 0;                       
wc = (2.*pi*c)/lambda;        
FWHM1 = 1e15*(lambda^2)/(c*bandw);  
offset = lambda - 2.*bandw;
start=0;
stop=4*bandw;
N=(stop-start)*fs;
w=linspace(start,stop,N);
w0=(stop-start)/2;
gau=exp(-4*log(2.)*((w-w0)/bandw).^2).*exp(1i*(wc*(w-w0)+ph));
ww = (w+offset)*1e9;                      
figure(1)
plot(ww,real(gau),ww,abs(gau))
xlabel('Wavelength (nm)')
ylabel('Power')
title('{\bf Periodogram}')
hold on
m = length(gau);          % Window
n = pow2(nextpow2(m));  % Transform
y = n*ifftshift(ifft(gau,n)); % DFT
t = (0:n-1)*(fs/n);     % Frequency
power = y.*conj(y)/n;   % Power 
figure(2)
plot(t,power,t,real(y))
xlabel('Time (s)')
ylabel('Intensity')
title('{\bf Time-Instensity Plot}')
hold on