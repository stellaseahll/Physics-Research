clear; clf;
TC = 1;
TH = 2;
TW = 8;
g = logspace(-3,0,50);
% % g(end) = [];
kappa = logspace(-4,-1,50);

Es = [5 1 4];
ks = [5 0 4];
% for i = 1:length(TW)
for j = 1:length(kappa)
    j
    for k = 1:length(g)
        
%             ks = ones(1,3)*kappa(j);
        M = model3QubitsFridgeOld(g(k),[TH TC TW],kappa(j)*ks,Es,'ohmic','x','p');
        M.findSS();
        QH(j,k) = M.getTotalHeatFlow(1);
        QC(j,k) = M.getTotalHeatFlow(2);
        QW(j,k) = M.getTotalHeatFlow(3);
        EH(j,k) = M.getEnergy(1);
        EC(j,k) = M.getEnergy(2);
        EW(j,k) = M.getEnergy(3);
%         C = M.getAllCoherence();
%         C1(j,k) = C(1);
%         C2(j,k) = C(2);
%         C3(j,k) = C(3);
%         E(j,k) = M.getEntanglement();
    end
end
% end
save nocoldbath.mat
