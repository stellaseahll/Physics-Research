TC = 1;
TH = 2;
TW = 8;
g = logspace(-4,0,100);
kappa = logspace(-4,-1,100);
Es = [5 1 4];
% QC = cell(1,length(TH));
% QH = cell(1,length(TH));
% QW = cell(1,length(TH));
% EC = cell(1,length(TH));
% EH = cell(1,length(TH));
% EW = cell(1,length(TH));
% for i = 1:length(TH)
%     i
    for j = 1:length(kappa)
       j
        for k = 1:length(g)
            ks = ones(1,3)*kappa(j).*Es;
            M = model3QubitsFridgeGM(g(k),[TH TC TW],ks,Es,'ohmic','x','p');
            M.findSS();
            QHtmp(j,k) = M.getTotalHeatFlow(1);
            QCtmp(j,k) = M.getTotalHeatFlow(2);
            QWtmp(j,k) = M.getTotalHeatFlow(3);
            EHtmp(j,k) = M.getEnergy(1);
            ECtmp(j,k) = M.getEnergy(2);
            EWtmp(j,k) = M.getEnergy(3);
        end
    end
    QH = QHtmp;
    QC = QCtmp;
    QW = QWtmp;
    EH = EHtmp;
    EC = ECtmp;
    EW = EWtmp;
% end
save ohmic_xxx_partial_TC1_TH2_TW8.mat
% clear;clc;
% 
% TC = 1;
% TH = linspace(1,5,100);
% TW = linspace(1,21,100);
% g = [0.1:0.1:0.9];
% Es = [5 1 4];
% ks = ones(1,3)*0.001.*Es;
% QC = cell(1,length(g));
% QH = cell(1,length(g));
% QW = cell(1,length(g));
% EC = cell(1,length(g));
% EH = cell(1,length(g));
% EW = cell(1,length(g));
% for i = 1:length(g)
%     tic
%     for j = 1:length(TH)
%         
%         for k = 1:length(TW)
%             M = model3QubitsFridge(g(i),[TH(j) TC TW(k)],ks,Es,'ohmic','x','p');
%             M.findSS();
%             QHtmp(j,k) = M.getTotalHeatFlow(1);
%             QCtmp(j,k) = M.getTotalHeatFlow(2);
%             QWtmp(j,k) = M.getTotalHeatFlow(3);
%             EHtmp(j,k) = M.getEnergy(1);
%             ECtmp(j,k) = M.getEnergy(2);
%             EWtmp(j,k) = M.getEnergy(3);
%         end
%     end
%     QH{i} = QHtmp;
%     QC{i} = QCtmp;
%     QW{i} = QWtmp;
%     EH{i} = EHtmp;
%     EC{i} = ECtmp;
%     EW{i} = EWtmp;
%     
%     fprintf('g %i -- ',i);
%     toc
% end
% save ohmic_xxx_partial_diffkappa_g(0.1,0.9).mat