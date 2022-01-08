%Compare quantum vs classical
Gamma = logspace(-2,2,21);
beta =  1;
dbeta = beta*0.001;
p0 = 1;
rho0 = diag([p0,1-p0]);
K = 0.1; %probability of NO swap
for i = 1:3
    i
    for j = 1:length(Gamma)
        j
        for k = 1:length(beta)
            M  = singlePassXY(K,Gamma(j),beta(k),dbeta(k),rho0,i);
            FQ{i}(j,k) = M.getFish();
            FC{i}(j,k) = getClassicalFish(K,Gamma(j),beta(k),dbeta(k),p0,i);
        end
    end
end
