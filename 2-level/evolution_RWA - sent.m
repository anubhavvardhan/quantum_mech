%Two-level Time evolution
function [r,rg] = evolution_RWA(Nf,gauf,mrabi,wc,hc,dt)
rho = zeros(2,2,Nf);
rho(:,:,1) = [[1,0];[0,0]];
rg = zeros(Nf,1);
r = zeros(Nf,1);
rg(1) = rho(1,1,1);
r(1) = rho(2,2,1);
%----Normalising gaussian to mrabi*pi area----
area=sum(abs(gauf)*dt);
E=mrabi*pi*gauf/area;
figure(4)
plot(E)
for iter=1:Nf-1
    H = 0.5*hc*[[-wc,E(iter)/hc];[conj(E(iter))/hc,wc]];
        [vec,val] = eig(H);
    rho(:,:,iter) = vec\rho(:,:,iter)*vec;
    e2 = [[exp(-1i*val(1,1)*dt), 0];[0,exp(-1i*val(2,2)*dt)]];
    e1 = [[exp(1i*val(1,1)*dt), 0];[0,exp(1i*val(2,2)*dt)]];
    rho(:,:,iter+1) = vec*e2*rho(:,:,iter)*e1/vec;
    r(iter+1) = rho(2,2,iter+1);
    rg(iter+1) = rho(1,1,iter+1);
end