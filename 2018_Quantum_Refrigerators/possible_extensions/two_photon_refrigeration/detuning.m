clc;clear;
TC = 1;
TH = 3;
TW = 5;
g = linspace(0,1,100);
kappa = [0.001 0.01];
delta = linspace(-1,1,100);
QC = cell(1,length(kappa));
QH = cell(1,length(kappa));
QW = cell(1,length(kappa));
EC = cell(1,length(kappa));
EH = cell(1,length(kappa));
EW = cell(1,length(kappa));
for i = 2:length(kappa)
    i
    for j = 1:length(delta)
       j
        for k = 1:length(g)
            Es = [5+delta(j) 1 4];
             
            M = model3QubitsFridge(g(k),[TH TC TW],ks,Es,'ohmic','x','p');
            M.findSS();
            QHtmp(j,k) = M.getTotalHeatFlow(1);
            QCtmp(j,k) = M.getTotalHeatFlow(2);
            QWtmp(j,k) = M.getTotalHeatFlow(3);
            EHtmp(j,k) = M.getEnergy(1);
            ECtmp(j,k) = M.getEnergy(2);
            EWtmp(j,k) = M.getEnergy(3);
        end
    end
    QH{i} = QHtmp;
    QC{i} = QCtmp;
    QW{i} = QWtmp;
    EH{i} = EHtmp;
    EC{i} = ECtmp;
    EW{i} = EWtmp;
end
save detuning_ohmic_xxx_partial.mat