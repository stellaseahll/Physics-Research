function [F rho rhoAnp]=getClassicalFish(gt,gammat,nbar,dnbar,pA,np,isSteady)
n = [nbar-dnbar, nbar,nbar+dnbar];
L = 1-exp(-gammat*(2*n+1));
K = (sin(gt))^2;
pGibbs = (n+1)./(2*n+1);
if isSteady
    pS = (L.*pGibbs + K*(1-L)*pA)./(1-(1-K)*(1-L));
else
    pS = pGibbs;
end
rhoA = [pA, 1-pA];
for i = 1:3
    rhoS{i} = [pS(i), 1-pS(i)];
    rho{i} = kron(rhoS{i},rhoA);
    for j = 2:np
        rho{i} = kron(rho{i},rhoA);
    end
end
for j = 1:np
    for i = 1:3
        rho{i} = processSwap(rho{i},j,np,K);
        rho{i} = processBath(rho{i},pGibbs(i),L(i));
    end
end

for i = 1:3
    rhoAnp{i} = rho{i}(1:(end/2)) + rho{i}((end/2+1):end);
end
drho = (rhoAnp{3}-rhoAnp{1})/2/dnbar;
F = sum((drho.^2)./rhoAnp{2});

end



function rho = processBath(rho,pGibbs,L)
p0 = rho(1:(end/2));
p1 = rho((end/2+1):end);
rho(1:(end/2)) = ((1-L)+pGibbs*L)*p0 + pGibbs*L*p1;
rho((end/2+1):end) = ((1-pGibbs)*L)*p0 + (1-L+(1-pGibbs)*L)*p1;
end

function rho = processSwap(rho,idx,np,K)
idx01 = find(kron(kron([1 0],kron(ones(1,2^(idx-1)),[0 1])),ones(1,2^(np-idx))));
idx10 = find(kron(kron([0 1],kron(ones(1,2^(idx-1)),[1 0])),ones(1,2^(np-idx))));
tmp1 = rho(idx01);
tmp2 = rho(idx10);
rho(idx01) = (1-K)*tmp1+K*tmp2;
rho(idx10) = (1-K)*tmp2+K*tmp1;
end