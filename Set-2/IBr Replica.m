clear
clc

hdash = 1.0*10^(-35);
%nu = 48.70388349;                          %for IBr
nu = 8.8698485d10;
di = 0.04;
popwt = 0.1;
tau = 20000;
ts = 500;
nt = tau/ts;

init_psi = dlmread('starting.txt');
nx = length(init_psi);

VG = dlmread('gdpotential.txt');
VE = dlmread('expotential.txt');
dr = VG(2,1) - VG(1,1);

rmin = 4;
rmax = 8;

E0 = 0.001;
wid_t = 1240.2412372725d0;                  %E = E0*exp(-x^2/wid_t)*exp(i*w0t)
w0 = 9.42637536d-3;

DG = 0.1025;
aG = 0.8172;
RG = 4.668;

DE = 0.03299;
aE = -0.7073;
RE = 5.166;

exgoal = sqrt(popwt) * init_psi(:,2);   % making exgoal 'exact replica' of ground state
prev_gket = init_psi(:,2);
present_gket = zeros(nx,1);
present_eket = zeros(nx,1);
prev_eket = zeros(nx,1);

prelaser = zeros(nt,1);
newlaser = zeros(nt,1);

for t = 1:nt
    overlap = 0;
    for n = 2:nx-1
        overlap = overlap + conj(prev_gket(n))*prev_eket(n)*dr;
        prelaser(t) = overlap;
        present_gket(n) = (prev_gket(n) + hdash*1i*ts/(nu*dr^2) * (prev_gket(n+1) - 2*prev_gket(n) + prev_gket(n-1)) - 2*1i*2*ts*di*E0*exp(-(4+n*dr)^2/wid_t)*cos(w0*t)*prev_eket(n)/hdash)/(1+2*ts*1i*DG*(exp(-2*aG*(n*dr+4-RG))-2*exp(-aG*(n*dr+4-RG)))/hdash);
        present_eket(n) = (prev_eket(n) + hdash*1i*ts/(nu*dr^2) * (prev_eket(n+1) - 2*prev_eket(n) + prev_eket(n-1)) - 2*1i*2*ts*di*E0*exp(-(4+n*dr)^2/wid_t)*cos(w0*t)*prev_gket(n)/hdash)/(1+2*ts*1i*DE*(exp(-2*aE*(n*dr+4-RE))-2*exp(-aE*(n*dr+4-RE)))/hdash);
    end
end


for p = 1:nt

    newlaser(p) = conj(prelaser(nt+1-p));
end

dlmwrite('laserpulsedata.txt',newlaser,'delimiter','\t','precision', 30);  % laserpulsedata will have the new laser  