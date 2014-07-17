function r = simulate()
w=2;
%wL= 9.42637536d-3;
wL=10;
wint = 10;
alpha=1;
%d = 0.04;
d=1;
pi=3.14159;
time = 2*pi;
nopoints = 5000;
Ef = zeros(nopoints);
t = linspace(0,time,nopoints);
dt = (time+1)/nopoints;
rho = zeros(2,2,nopoints);
r = zeros(nopoints,1);
rho(:,:,1) = [[1,0];[0,0]];
r(1) = rho(2,2,1);
for iter = 1:(nopoints-1)
    [H,E] = hamiltonian(t(iter),time/2,w,wL,d,alpha,wint);
    Ef(iter) = E;
    [vec,val] = eig(H);
    e2 = [[exp(-1i*val(1,1)*dt), 0];[0,exp(-1i*val(2,2)*dt)]];
    e1 = [[exp(1i*val(1,1)*dt), 0];[0,exp(1i*val(2,2)*dt)]];
    rho(:,:,iter+1) = vec*e2*vec\rho(:,:,iter)*vec*e1/vec;
    r(iter+1) = rho(2,2,iter+1);
end
figure(1)
plot(t, Ef)
figure(2)
plot(t,abs(r))