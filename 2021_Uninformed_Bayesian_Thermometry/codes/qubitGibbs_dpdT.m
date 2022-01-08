function p = qubitGibbs_dpdT(T)
%helper function to compute derivative of excitation probability of qubit w.r.t. T, 
%given T in units of hbar*omega/kB

p = 1./(exp(1./T) + exp(-1./T) + 2)./T.^2;

p(T==0) = 0;
p(isinf(T)) = 0;