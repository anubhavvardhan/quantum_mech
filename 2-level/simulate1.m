function x=simulate()
w=2;
%wL= 9.42637536d-3;
wL=10;
wint = 10;
alpha=1;
%d = 0.04;
d=1;
pi=3.14159;
time = 2*pi;
nopoints = 3;
t = linspace(0,time,nopoints);
dt = (time+1)/nopoints;
rho = zeros(2,2,nopoints);
rho(:,:,1) = [[1,0];[0,0]];
for iter = 1:(nopoints-1)
    H = hamiltonian(t(iter),time/2,w,wL,d,alpha,wint);
    [vec,val] = eig(H);
    e2 = [[exp(-i*val(1,1)*dt), 0];[0,exp(-i*val(2,2)*dt)]];
    e1 = [[exp(i*val(1,1)*dt), 0];[0,exp(i*val(2,2)*dt)]];
    rho(:,:,iter+1) = vec*e2*inv(vec)*rho(:,:,iter)*vec*e1*inv(vec);
end