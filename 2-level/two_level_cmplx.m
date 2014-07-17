clear all
[lambda bandw c hbar ph wc gaun mrabi FWHM1 gauf gauf1 dt ww Nf tt] = variable;
% mrabi=50.;
%----------------------------------------------------------------------------------
figure(1)
plot(ww, real(gaun), ww, imag(gaun), ww, abs(gaun));
xlabel('Wavelength (nm)');
ylabel('Intensity (arbitrary units)');
%----------------------------------------------------------------------------------
figure(2)
plot(tt*1e15, real(gauf), tt*1e15, abs(gauf)), grid on;
xlabel('Time (fs)');
ylabel('Intensity (arbitrary units)');
title(sprintf('Pulse in Time Domain with FWHM = %d fs',FWHM1));
%----------------------------------------------------------------------------------
nct = round((FWHM1*1e-15)*wc/(2.*pi));     % Number of cycles in time FWHM
display(FWHM1)
display(nct)
%----------------------------------------------------------------------------------
% [r,rg] = evolution(Nf,gauf1,mrabi,wc,hbar,dt);
[r,rg] = evolution(Nf,gauf,mrabi,wc,hbar,dt);
figure(4)
plot(tt*1e15,real(r),tt*1e15,abs(rg),tt*1e15,abs(gauf))
xlabel('Time (fs)');
ylabel('Population');
legend excited ground pulse
title(sprintf('Population Evolution in Time Domain with %d pi Pulse Area',mrabi));
%----------------------------------------------------------------------------------