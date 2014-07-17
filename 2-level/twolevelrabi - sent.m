%---------Two Level System---------------
%-------setting gaussian pulse-----------
FWHM=0.5;
ome = 100;
ph = 0;
delta = 0.;
start=0;
stop=4*FWHM;
dt=0.01;
N=(stop-start)/dt;
t=linspace(start,stop,N);
t0=(stop-start)/2;
field_p = exp(1i*(ome*t + ph));
gaun=exp(-4*log(2)/(FWHM^2)*(t-t0).^2);
%-------Rabi Frequency Step---------------
mrabi=4;
dOm=.05;
N_rabi=mrabi/dOm;
Om=linspace(0,mrabi,N_rabi);
%----Normalising gaussian to m*pi area----
m=1;
area=sum(gaun)*dt;
gaun=m*pi*gaun.*field_p/area;
%-----------------------------------------
figure(1);
plot(t, abs(gaun),t, real(gaun))
xlabel('Time(ps)');
ylabel('Intensity(arbitrary units)');
%-----------------------------------------
psi1=zeros(N_rabi);
psi2=zeros(N_rabi);    

for r=1:1:N_rabi
    ro = [1 0;
          0 0];
for j = 1:1:N
    %Hamiltonian
    H = [delta Om(r)*gaun(j)/2;
        Om(r)*gaun(j)/2 0]; 
    
    [VH,DH] = eig(H);
    
    ro=VH\ro*VH;
    u1=expm(1i*DH*dt);
    u2=expm(-1i*DH*dt);
    ro=u2*ro*u1;
    ro=VH*ro/VH;
    
end
psi1(r)=abs(ro(1,1));
psi2(r)=abs(ro(2,2));
end
figure(2);
plot(Om,psi1,'b',Om,psi2,'g');
xlabel('Rabi Frequency(pi units)');
ylabel('Intensity(arbitrary units)');

%------------------------------------------

  





