% clear;clc;
% k = 0.001;
% g = 0.001;
% EH = 4.5:0.5:10;
% EC = 1;
% EW = EH-EC;
% TW = 4.5:0.5:20;
% TH = 2;
% TC = 1;
% r = 0;
% for i = 1:length(EH)
%     for j = 1:length(TW)
%         E = [EH(i) EC EW(i)];
%         T = [TH TC TW(j)];
%         % One photon
%         M = model3SpinsFridge(g,T,k*E,E,'f','r','l',1);
%         M.findSS();
%         qc = M.getTotalHeatFlow(2);
%         % Two photons (WC)
%         nh = 1./(exp(E(1)./T(1))-1);
%         nc = 1./(exp(2*E(2)./T(2))-1);
%         nw = 1./(exp(2*E(3)./T(3))-1);
% 
%         rhoh = exp(-E(1)./T(1)*(0:2));
%         rhoh = diag(rhoh/sum(rhoh));
%         N = model3SpinsFridgeTwoPhotonWC(g,T,k*E,E,'f','r','l',1);
%         rhow = diag([1 0 0]);
%         rhoc = diag([0 1 0]);
%         rho2 = kron(rhoc,rhow);    
%         rhocw = rho2;
%         rhocw = rhocw/trace(rhocw);
%         N.findSS(kron(rhoh,rhocw));    
%         QC(i,j) = N.getTotalHeatFlow(2);
%         QCrat(i,j) = (QC(i,j)-qc)/qc;
%     end
% end
% clear M; clear N;
% save getMaxImprovement.mat;
% load save getMaxImprovement.mat;
% contourf(TW,EH,QCrat)

