% Calculate Probability Density as a function of time
% for a particle trapped in a double-well potential
%
% Potential due to two square wells of width 2w
% and a distance 2a apart
w = L/50; a = 3*w;
U = -100*( heaviside(x+w-a) - heaviside(x-w-a) ...
+ heaviside(x+w+a) - heaviside(x-w+a));
% Finite-difference representation of Laplacian and Hamiltonian,
% where hbar = m = 1.
e = ones(N,1); Lap = spdiags([e -2*e e],[-1 0 1],N,N)/dx^2;
H = -(1/2)*Lap + spdiags(U,0,N,N);
% Find and sort lowest nmodes eigenvectors and eigenvalues of H
nmodes = 2; options.disp = 0; [V,E] = eigs(H,nmodes,’sa’,options);
[E,ind] = sort(diag(E));% convert E to vector and sort low to high
V = V(:,ind); % rearrange coresponding eigenvectors
% Rescale eigenvectors so that they are always
% positive at the center of the right well
for c = 1:nmodes
V(:,c) = V(:,c)/sign(V((3*N/4),c));
end
%****************************************************************
% Compute and display normalized prob. density rho(x,T)
%****************************************************************
% Parameters for solving the problem in the interval 0 < T < TF
TF = 10; % Length of time interval
NT = 100; % No. of time points
T = linspace(0,TF,NT); % Time vector
% Compute probability normalisation constant (at T=0)
psi_o = 0.5*V(:,1)+0.5*V(:,2); % wavefunction at T=0
sq_norm = psi_o’*psi_o*dx; % square norm = |<ff|ff>|^2
Usc = U*max(abs(V(:)))/max(abs(U)); % rescale U for plotting
% Compute and display rho(x,T) for each time T
for t=1:NT; % time index parameter for stepping through loop
% Compute wavefunction psi(x,T) and rho(x,T) at T=T(t)
psi = 0.5*V(:,1)*exp(-1i*E(1)*T(t)) ...
+ 0.5*V(:,2)*exp(-1i*E(2)*T(t));
rho = conj(psi).*psi/sq_norm; % normalized probability density
% Plot rho(x,T) and rescaled potential energy Usc
plot(x,rho,’o-k’, x, Usc,’.-b’); axis([-L/8 L/8 -1 6]);
lgnd_str = [repmat(’T = ’,1,1),num2str(T)];
text(-0.12,5.5,lgnd_str, ’FontSize’, 18);
xlabel(’x [m]’, ’FontSize’, 24);
ylabel(’probability density [1/m]’,’FontSize’, 24);
pause(0.05); % wait 0.05 seconds
end