function rho = prepareTwoQubitState(r,phi,p1,t1,p2,t2)

ket1p = [cos(t1/2); sin(t1/2)*exp(1i*p1)];
ket1m = [sin(t1/2)*exp(-1i*p1); -cos(t1/2)];
ket2p = [cos(t2/2); sin(t2/2)*exp(1i*p2)];
ket2m = [sin(t2/2)*exp(-1i*p2); -cos(t2/2)];

ket = sqrt(r) * kron(ket1p,ket2p) + exp(1i*phi)*sqrt(1-r) * kron(ket1m,ket2m);
rho = ket*ket';