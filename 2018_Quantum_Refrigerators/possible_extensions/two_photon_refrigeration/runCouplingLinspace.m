clear; 
TC = 1;
TH = 1.2;
TW = 4;
Es = [3 1 2];
fnam = 'runLog200_O-X-P_E3-1-2_T1.2-1-4.mat';

N = 200;

%g = (1:N)/N*1.0;
%kappa = (1:N)/N*0.1;
g = logspace(-3,0,N+1);
kappa = logspace(-4,-1,N+1);


% for i = 1:length(TW)
for j = 1:length(kappa)
    tic
    fprintf('chi step %i -- ',j)
    for k = 1:length(g)
        ks = kappa(j).*Es;
%             ks = ones(1,3)*kappa(j);
        M = model3QubitsFridgeOld(g(k),[TH TC TW],ks,Es,'ohmic','x','p');
        M.findSS();
        rhoD = eig(M.rho);
        minEig(j,k) = min(rhoD);
        trRho(j,k) = sum(rhoD);
        QH(j,k) = M.getTotalHeatFlow(1);
        QC(j,k) = M.getTotalHeatFlow(2);
        QW(j,k) = M.getTotalHeatFlow(3);
        EH(j,k) = M.getEnergy(1);
        EC(j,k) = M.getEnergy(2);
        EW(j,k) = M.getEnergy(3);
        C = M.getAllCoherence();
        C1(j,k) = C(1);
        C2(j,k) = C(2);
        C3(j,k) = C(3);
        p100(j,k) = M.rho(5,5);
        p011(j,k) = M.rho(4,4);
        S(j,k) = M.getEntropy(M.rho);
        SH(j,k) = M.getEntropy(M.rhoH);
        SC(j,k) = M.getEntropy(M.rhoC);
        SW(j,k) = M.getEntropy(M.rhoW);
        SCW(j,k) = M.getEntropy(M.rhoCW);
        SHW(j,k) = M.getEntropy(M.rhoHW);
        SHC(j,k) = M.getEntropy(M.rhoHC);
        negH(j,k) = M.getNegativityHot();
        negC(j,k) = M.getNegativityCold();
    end
    toc
end
% end
save(fnam)


%% Check parameters Brunner2014

clear; clf;
TC = 1;
TH = 1.1;
TW = 1e4;
Es = [302 2 300];

N = 400;

%g = (1:N)/N*1.0;
%kappa = (1:N)/N*0.1;
g = logspace(-4,0,N+1);
kappa = logspace(-5,-1,N+1);

% for i = 1:length(TW)
for j = 1:length(kappa)
    tic
    fprintf('chi step %i -- ',j)
    for k = 1:length(g)
        ks = kappa(j).*Es;
%             ks = ones(1,3)*kappa(j);
        M = model3QubitsFridgeOld(g(k),[TH TC TW],ks,Es,'ohmic','x','p');
        M.findSS();
        QH(j,k) = M.getTotalHeatFlow(1);
        QC(j,k) = M.getTotalHeatFlow(2);
        QW(j,k) = M.getTotalHeatFlow(3);
        EH(j,k) = M.getEnergy(1);
        EC(j,k) = M.getEnergy(2);
        EW(j,k) = M.getEnergy(3);
        C = M.getAllCoherence();
        C1(j,k) = C(1);
        C2(j,k) = C(2);
        C3(j,k) = C(3);
        p100(j,k) = M.rho(5,5);
        p011(j,k) = M.rho(4,4);
        S(j,k) = M.getEntropy(M.rho);
        SH(j,k) = M.getEntropy(M.rhoH);
        SC(j,k) = M.getEntropy(M.rhoC);
        SW(j,k) = M.getEntropy(M.rhoW);
        SCW(j,k) = M.getEntropy(M.rhoCW);
        SHW(j,k) = M.getEntropy(M.rhoHW);
        SHC(j,k) = M.getEntropy(M.rhoHC);
        negH(j,k) = M.getNegativityHot();
        negC(j,k) = M.getNegativityCold();
        %resonant Hamiltonian
        M = model3QubitsFridgeOld(g(k),[TH TC TW],ks,Es,'ohmic','r','p');
        M.findSS();
        QHr(j,k) = M.getTotalHeatFlow(1);
        QCr(j,k) = M.getTotalHeatFlow(2);
        QWr(j,k) = M.getTotalHeatFlow(3);
        EHr(j,k) = M.getEnergy(1);
        ECr(j,k) = M.getEnergy(2);
        EWr(j,k) = M.getEnergy(3);
        C = M.getAllCoherence();
        C1r(j,k) = C(1);
        C2r(j,k) = C(2);
        C3r(j,k) = C(3);
        p100r(j,k) = M.rho(5,5);
        p011r(j,k) = M.rho(4,4);
        Sr(j,k) = M.getEntropy(M.rho);
        SHr(j,k) = M.getEntropy(M.rhoH);
        SCr(j,k) = M.getEntropy(M.rhoC);
        SWr(j,k) = M.getEntropy(M.rhoW);
        SCWr(j,k) = M.getEntropy(M.rhoCW);
        SHWr(j,k) = M.getEntropy(M.rhoHW);
        SHCr(j,k) = M.getEntropy(M.rhoHC);
        negHr(j,k) = M.getNegativityHot();
        negCr(j,k) = M.getNegativityCold();
        %local ME
        M = model3QubitsFridgeOld(g(k),[TH TC TW],ks,Es,'ohmic','x','l');
        M.findSS();
        QHl(j,k) = M.getTotalHeatFlow(1);
        QCl(j,k) = M.getTotalHeatFlow(2);
        QWl(j,k) = M.getTotalHeatFlow(3);
        EHl(j,k) = M.getEnergy(1);
        ECl(j,k) = M.getEnergy(2);
        EWl(j,k) = M.getEnergy(3);
        C = M.getAllCoherence();
        C1l(j,k) = C(1);
        C2l(j,k) = C(2);
        C3l(j,k) = C(3);
        p100l(j,k) = M.rho(5,5);
        p011l(j,k) = M.rho(4,4);
        Sl(j,k) = M.getEntropy(M.rho);
        SHl(j,k) = M.getEntropy(M.rhoH);
        SCl(j,k) = M.getEntropy(M.rhoC);
        SWl(j,k) = M.getEntropy(M.rhoW);
        SCWl(j,k) = M.getEntropy(M.rhoCW);
        SHWl(j,k) = M.getEntropy(M.rhoHW);
        SHCl(j,k) = M.getEntropy(M.rhoHC);
        negHl(j,k) = M.getNegativityHot();
        negCl(j,k) = M.getNegativityCold();
        %resonant Hamiltonian & local model
        M = model3QubitsFridgeOld(g(k),[TH TC TW],ks,Es,'ohmic','r','l');
        M.findSS();
        QHrl(j,k) = M.getTotalHeatFlow(1);
        QCrl(j,k) = M.getTotalHeatFlow(2);
        QWrl(j,k) = M.getTotalHeatFlow(3);
        EHrl(j,k) = M.getEnergy(1);
        ECrl(j,k) = M.getEnergy(2);
        EWrl(j,k) = M.getEnergy(3);
        C = M.getAllCoherence();
        C1rl(j,k) = C(1);
        C2rl(j,k) = C(2);
        C3rl(j,k) = C(3);
        p100rl(j,k) = M.rho(5,5);
        p011rl(j,k) = M.rho(4,4);
        Srl(j,k) = M.getEntropy(M.rho);
        SHrl(j,k) = M.getEntropy(M.rhoH);
        SCrl(j,k) = M.getEntropy(M.rhoC);
        SWrl(j,k) = M.getEntropy(M.rhoW);
        SCWrl(j,k) = M.getEntropy(M.rhoCW);
        SHWrl(j,k) = M.getEntropy(M.rhoHW);
        SHCrl(j,k) = M.getEntropy(M.rhoHC);
        negHrl(j,k) = M.getNegativityHot();
        negCrl(j,k) = M.getNegativityCold();
    end
    toc
end
% end
save runBrunner400partialX.mat