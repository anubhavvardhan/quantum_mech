% Symmetric Potential
% define matrix that is discretized
% version of Laplacian, ie, del squared operator
% ----------------------------------------------
n=100; deltax=1/n; x=0:deltax:1;
delsq = toeplitz([-2 1 zeros(1,n-1)]);
%not including 1/deltax^2 term
% define potential function
% -------------------------
Vmax=.1; v1=0*ones(1,40); v2=Vmax*sin((0:20)/20*pi).^2;
v3=0*ones(1,40); v=[v1,v2,v3];
do=1; 
if do
% hamiltonian operator
% --------------------
hbar=1; m=1; H = -hbar^2/(2*m)*delsq + diag(v); A = 1/(1i*hbar)*H;
% finding energies of eigenstates
[V,D]=eig(H);
[E,index]=sort(diag(D)); %sort energies
V=V(:,index); %sort eigenfunctions
for m=1:(n+1); %normalizing eigenfunctions
V(:,m)=V(:,m)/norm(V(:,m));
end;
% initial wave fct
% ----------------
% stationary wave packet;
% Psi0 = exp(-(x-0.2).^2/0.1^2);
% moving wave packet;
k=2*pi/0.2;
Psi0=exp(-(x-0.2).^2/0.1^2).*exp(1i*k*x); Psi0=conj(Psi0');
Psi0 = Psi0/norm(Psi0); meanE=Psi0'*H*Psi0; Alpha=V'*Psi0;
absAlpha=abs(Alpha).^2;
fprintf('Mean energy is %f\n', meanE);
fprintf('Vmax is %f \n\n', Vmax);
% time evolution
% --------------
deltat=40; M=expm(deltat*A); PsiMat=Psi0; Psi=Psi0;
nsteps=500; for z=1:nsteps; Psi=M*Psi; PsiMat=[PsiMat Psi]; end;
PDense=abs(PsiMat).^2; % Probability Density.
PDensel=PDense(1:50,:); Pl=sum(PDensel);
%deltat=200; M=expm(deltat*A);
%PsiMat=Psi0; Psi=Psi0;
%nsteps=10;
%for z=1:nsteps; Psi=M*Psi; PsiMat=[PsiMat Psi]; end;
%PDense2=abs(PsiMat).^2; % Probability Density.
end;
% some plots
% ----------
close all; figure(1); plot(x,v);
xlabel('x'); ylabel('V'); title('potential function');
grid on;
if do
figure(2);
m=4;
for l=1:m
subplot(m,1,l);
plot(x, real(V(:,l)),x,v/Vmax*.3,':');
axis([0 1 -.3 .3]);
end;
subplot(m,1,1); title(' Low energy eigenfunctions');
subplot(m,1,m); xlabel('x');
figure(3);
subplot(2,1,1);
plot(x,PDense(:,1)); grid on;
xlabel('x'); ylabel('|Psi0|^2');
title('Wavefunction at t=0');
subplot(2,1,2);
plot(0:100,absAlpha,'o',0:100,absAlpha);
grid on;
xlabel('n');
ylabel('|Alphasubn|^2');
title('|Alphasubn|^2 vs. n');
figure(4);
m=8;
for l=1:m;
subplot(m,1,l);
plot(x,PDense(:,l),x,v,':');
axis off;
end;
subplot(m,1,1); title('Psi at intervals of 40 t.u.');
subplot(m,1,m); xlabel('x');
figure(5);
m=8;
for l=1:m;
subplot(m,1,l);
plot(x,PDense(:,l+8),x,v,':');
axis off;
end;
subplot(m,1,1); title('Psi at intervals of 40 t.u.');
subplot(m,1,m); xlabel('x');
figure(6); plot(0:nsteps,Pl)
grid on; zoom on;
axis([0 nsteps 0 1]);
xlabel('time'); ylabel('Pleft');
title('Pleft vs. t');
% momentum operator
% -----------------
%P= hbar/(2*deltax)*toeplitz([0, i,zeros(1,n-1)]);
end;