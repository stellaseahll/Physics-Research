p = 0.99:-0.01:0.01;
H = [0;1];
beta = 1;
alpha = 2.8;
for i = 1:length(p)
    rho(:,i) = fsolve(@(rho) getRho(rho,[p(i);1-p(i)],alpha,beta,H), [p(i);1-p(i)]);
end
