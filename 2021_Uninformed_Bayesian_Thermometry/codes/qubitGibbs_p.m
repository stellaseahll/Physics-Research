function p = qubitGibbs_p(T)
%helper function to compute excitation probability of qubit, given T in
%units of hbar*omega/kB

p = 1./(exp(1./T)+1);
p(T==0) = 0;
p(isinf(T)) = 0.5;