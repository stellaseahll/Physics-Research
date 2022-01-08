function [rhoP,rhoM] = getProjectedState(rho)
dp = length(rho);
[Pp Pm] = getIntegral(dp);
rhoP = zeros(dp);
rhoM = zeros(dp);
for j = 1:dp
    for k = 1:dp
        for m = 1:dp
            for n = 1:dp
                rhoP(j,k) = rhoP(j,k)+rho(m,n)*Pp(m,j)*Pp(n,k);
                rhoM(j,k) = rhoM(j,k)+rho(m,n)*Pm(m,j)*Pm(n,k);
            end
        end
    end
end
    

end