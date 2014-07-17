function H = hamiltonian(t,centre,w,wL,d,alpha,wint)
E = exp(-alpha*(t-centre)*(t-centre))*cos(wint*(t-centre));
H = [[-w/2,-2*d*E*cos(wL*t)];[-2*d*E*cos(wL*t),w/2]];
