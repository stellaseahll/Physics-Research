function C = getRelCoh(rho)
%compute relative entropy of coherence
C = getSvn(diag(diag(rho))) - getSvn(rho);
end

function S = getSvn(rho)
k = eig(rho);
S = -k.*log(k);
S = sum(S(k>eps));
end