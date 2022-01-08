function d = getQ(rho,L)
    d = L*rho*L' - 0.5*(rho*L'*L + L'*L*rho);
end