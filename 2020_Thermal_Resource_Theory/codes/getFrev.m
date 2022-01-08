function Frev = getFrev(rhoi,rhof,Hs,beta)
%This function computes the change in free energy by a reversible process

%% First Step: Change Hamiltonian such that rhoi is a thermal state with respect to beta
H1 = -log(rhoi)/beta;
dE1 = sum(rhoi.*(H1-Hs));

%% Second Step: Change Hamiltonian adiabatically from H1->H2 such that we get a thermal state rhof at the end
H2 = -log(rhof)/beta;
dE2 = sum(rhof.*H2 - rhoi.*H1);
dS2 = -sum(rhof.*log(rhof)-rhoi.*log(rhoi));

%% Third Step: Change Hamiltonian back to original Hamiltonian
dE3 = sum(rhof.*(Hs-H2));

%% Net
Frev = dE1+dE2+dE3-dS2/beta;

end

