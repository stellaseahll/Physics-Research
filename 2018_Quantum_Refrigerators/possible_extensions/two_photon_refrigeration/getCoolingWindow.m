% clear;clc;
% k = 0.001;
% g = 0.001;
% TC = 1;
% TH = linspace(1,5,50);
% TW = linspace(1,21,50);
% E = [5 1 4];
% r1 = [1:-0.1:0; 0:0.1:1];
% QC = cell(1,length(r1));
% QH = cell(1,length(r1));
% QW = cell(1,length(r1));
% EC = cell(1,length(r1));
% EH = cell(1,length(r1));
% EW = cell(1,length(r1));

for i = 2:length(r1)
    i
    for j = 1:length(TH)
        j
        for k = 1:length(TW)
            k
            nh = 1./(exp(E(1)./TH(j))-1);
            nc = 1./(exp(2*E(2)./TC)-1);
            nw = 1./(exp(2*E(3)./TW(k))-1);
            rhow = [r1(1,i)*(nw+1)/(2*nw+1) 1-r1(1,i) r1(1,i)*(nw)/(2*nw+1)];
            rhoc = [r1(2,i)*(nc+1)/(2*nc+1) 1-r1(2,i) r1(2,i)*(nc)/(2*nc+1)];
            rhoh = exp(-E(1)./TH(j)*(0:2));
            rhoh = diag(rhoh/sum(rhoh));    
            rhow = diag(rhow);
            rhoc = diag(rhoc);
            N = model3SpinsFridgeTwoPhotonWC(g,[TH(j) TC TW(k)],k*E,E,'f','r','l',1);
            N.findSS(kron(rhoh,kron(rhoc,rhow))); 
            QHtmp(j,k) = N.getTotalHeatFlow(1);
            QCtmp(j,k) = N.getTotalHeatFlow(2);
            QWtmp(j,k) = N.getTotalHeatFlow(3);
%             EHtmp(j,k) = N.getEnergy(1);
%             ECtmp(j,k) = N.getEnergy(2);
%             EWtmp(j,k) = N.getEnergy(3);
        end
    end
    QH{i} = QHtmp;
    QC{i} = QCtmp;
    QW{i} = QWtmp;
%     EH{i} = EHtmp;
%     EC{i} = ECtmp;
%     EW{i} = EWtmp;
end

save coolingwindow_TC1_EH5_EW4_EC1.mat
