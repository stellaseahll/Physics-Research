function T = qubitGibbs_dTdp(p)
%helper function to compute derivative of T of qubit w.r.t. p, given excitation prob p in 0..1, 
%T in units of hbar*omega/kB

%T = 1./log(1./p - 1);
T = 1./log(1./p-1).^2 ./(p.*(1-p));

T(p==0) = inf;
T(p==1) = inf;
T(p==0.5) = inf;