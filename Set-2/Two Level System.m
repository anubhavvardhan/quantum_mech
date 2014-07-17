clear all
[lambda bandw c hbar ph wc gaun mrabi FWHM1 gauf gauf1 dt ww Nf tt] = variable;
%
figure(1)
plot(ww, real(gaun), ww, imag(gaun), ww, abs(gaun));
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
figure(2)
plot(tt*1e15, real(gauf), tt*1e15, abs(gauf)), grid on;
xlabel('Time (fs)');
ylabel('Intensity (arbitrary units)');
title(sprintf('Pulse in Time Domain with FWHM = %d fs',FWHM1));
cycles = round((FWHM1*1e-15)*wc/(2.*pi));     
display(FWHM1)
display(cycles)
[r,rg] = Evolution(Nf,gauf,mrabi,wc,hbar,dt);
figure(4)
plot(tt*1e15,real(r),tt*1e15,abs(rg),tt*1e15,abs(gauf))
xlabel('Time (fs)');
ylabel('Population');
legend excited ground pulse
title(sprintf('Population Evolution in Time Domain with %d pi Pulse Area',mrabi));