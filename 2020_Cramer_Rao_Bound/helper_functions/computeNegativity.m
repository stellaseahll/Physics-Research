function f = computeNegativity(rho)

rho(1:2,1:2) = rho(1:2,1:2).';
rho(3:4,3:4) = rho(3:4,3:4).';
rho(1:2,3:4) = rho(1:2,3:4).';
rho(3:4,1:2) = rho(3:4,1:2).';

R = eig(rho);
f = -sum(R(R<0));