% clear;clc;
% k = 0.001;
% g = 0.001;
% E = [5 1 4];
% TW = 3:30;
% TH = 2;
% TC = 1;
% r = 0:0.01:1;
% for i = 1:length(TW)
%     T = [TH TC TW(i)];
%     % One photon
%     M = model3SpinsFridge(g,T,k*E,E,'f','r','l',1);
%     M.findSS();
%     qc = M.getTotalHeatFlow(2);
% 
%     % Two photons (WC)
%     nh = 1./(exp(E(1)./T(1))-1);
%     nc = 1./(exp(2*E(2)./T(2))-1);
%     nw = 1./(exp(2*E(3)./T(3))-1);
% 
%     rhoh = exp(-E(1)./T(1)*(0:2));
%     rhoh = diag(rhoh/sum(rhoh));
% 
% 
%     N = model3SpinsFridgeTwoPhotonWC(g,T,k*E,E,'f','r','l',1);
% 
%     rhow = diag([1 0 0]);
%     rhoc = diag([1 0 0]);
%     rho1 = kron(rhoc,rhow);
% 
%     rhow = diag([1 0 0]);
%     rhoc = diag([0 1 0]);
%     rho2 = kron(rhoc,rhow);
%     for j = 1:length(r)
%         rhocw = r(j)*rho1 + (1-r(j))*rho2;
%         rhocw = rhocw/trace(rhocw);
%         N.findSS(kron(rhoh,rhocw));    
%         QC(i,j) = N.getTotalHeatFlow(2);
%     end
%     QCrat(i,:) = (QC(i,:)-qc)/qc;
% end
% save linearCombineSubspace_VaryTW_TC1_TH2_EH5_EW4_EC1.mat;
load linearCombineSubspace_VaryTW_TC1_TH2_EH5_EW4_EC1.mat;
hold on;
figure;
plot(r,QCrat);
figure;
plot(TW,QCrat(:,1))
