function f = getRho(rho,tau,alpha,beta,H)

lntau = log(tau);
lntau(tau<=0) = 0;
k = (1-(1-alpha)*(-lntau-log(sum(exp(-beta*H)))-beta*sum(rho.^alpha.*(H))/sum(rho.^alpha))).^(1/(1-alpha));
k = real(k);
k = k/sum(k);
f = k-rho;