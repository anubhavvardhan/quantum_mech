dat = dlmread('D:/Lab/ML null_re');          %data file
raw = 0.00003*dat;
autocorr = xcorr(raw(:,1));
Fs = 10^3;                                                                   %sampling frequency
NFFT = 2^nextpow2(length(autocorr));
Y = fft(autocorr(:),NFFT)/length(autocorr);
f = Fs/2*linspace(0,1,NFFT/2);
y = log10(abs(Y(1:NFFT/2)));
x = log10(f);
%plot(x,y);
yy = smooth(y,99);
plot(x,yy);
P = [x;yy']';
dlmwrite('D:/Lab/MLsmooth',P);