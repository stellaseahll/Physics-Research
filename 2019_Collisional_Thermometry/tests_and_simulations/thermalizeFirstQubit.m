function rho = thermalizeFirstQubit(rho0,gamma_t,nbar)
%Function applies local thermalizing channel to system qubit, assuming the
%total state rho0 is a (2*dP * 2*dP) matrix constructed using system qubit
%as first tensor factor.
% gamma_t: thermalization time gamma*t
% nbar: average excitation of bath
%Output is state matrix rho, works also with sparse matrices

%total probe dimension
dP = length(rho0)/2;

E = exp(-gamma_t*(2*nbar+1)/2);
E2 = exp(-gamma_t*(2*nbar+1));
p = (nbar+1)/(2*nbar+1);

%decay of coherences
rho = E*rho0;

%mixing of excited and ground populations
rho(1:dP,1:dP) = (p+(1-p)*E2) * rho0(1:dP,1:dP) + p*(1-E2) * rho0(dP+1:end,dP+1:end);
rho(dP+1:end,dP+1:end) = (1-p)*(1-E2) * rho0(1:dP,1:dP) + (1-p+p*E2) * rho0(dP+1:end,dP+1:end);