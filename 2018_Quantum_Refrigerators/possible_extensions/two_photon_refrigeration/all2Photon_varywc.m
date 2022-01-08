k = 0.001;
g = 0.01;
Ec = [linspace(1e-3,1.8374,245) logspace(log10(1.8374),log10(1.8746),100) linspace(1.8747,1.875,100)];
TW = 5;
TH = 2;
TC = 1;
T = [TH TC TW];
rho0 = zeros(3);
rho0(1,1) = 1;
rho1 = zeros(3);
rho1(2,2) = 1;
rho2 = zeros(3);
rho2(3,3) = 1;
rho{1} = kron(rho0,kron(rho0,rho0));
rho{2} = kron(rho0,kron(rho0,rho1));
rho{3} = kron(rho0,kron(rho1,rho0));
rho{4} = kron(rho0,kron(rho1,rho1));
for i = 245:length(Ec)
    i
    E = [5,Ec(i),5-Ec(i)];
    Nl = model3SpinsFridgeTwoPhoton(g,T,k*E,E,'o','r','l',1);
    Ng = model3SpinsFridgeTwoPhoton(g,T,k*E,E,'o','r','g',1);
    Np = model3SpinsFridgeTwoPhoton(g,T,k*E,E,'o','r','p',1);
    for j = 1:4
        Nl.findSS(rho{j});   
        Ng.findSS(rho{j});   
        Np.findSS(rho{j});
        QHl(i,j) = Nl.getTotalHeatFlow(1);
        QHg(i,j) = Ng.getTotalHeatFlow(1);
        QHp(i,j) = Np.getTotalHeatFlow(1);
        QCl(i,j) = Nl.getTotalHeatFlow(2);
        QCg(i,j) = Ng.getTotalHeatFlow(2);
        QCp(i,j) = Np.getTotalHeatFlow(2);
        QWl(i,j) = Nl.getTotalHeatFlow(3);
        QWg(i,j) = Ng.getTotalHeatFlow(3);
        QWp(i,j) = Np.getTotalHeatFlow(3);
    end
    save 2photonHWC_lgp_Varywc_g1e-2_TC1_TH2_TW5.mat;
%     QCrat(i,:) = (QC(i,:)-qc)/qc;
end