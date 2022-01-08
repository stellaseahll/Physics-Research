clear;clc;
k = 0.001;
g = 0.001;
E = [5 1 4];
T = [2 1 8];

% One photon
M = model3SpinsFridge(g,T,k*E,E,'f','r','l',1);
M.findSS();
qc = M.getTotalHeatFlow(2);

% Two photons (WC)
nh = 1./(exp(E(1)./T(1))-1);
nc = 1./(exp(2*E(2)./T(2))-1);
nw = 1./(exp(2*E(3)./T(3))-1);

rhoh = exp(-E(1)./T(1)*(0:2));
rhoh = diag(rhoh/sum(rhoh));

Nt = 1000;
r1 = rand(2,Nt);
N = model3SpinsFridgeTwoPhotonWC(g,T,k*E,E,'f','r','l',1);
for j = 1:Nt
    rhow = [r1(1,j)*(nw+1)/(2*nw+1) 1-r1(1,j) r1(1,j)*(nw)/(2*nw+1)];
    rhoc = [r1(2,j)*(nc+1)/(2*nc+1) 1-r1(2,j) r1(2,j)*(nc)/(2*nc+1)];
    rhow = diag(rhow);
    rhoc = diag(rhoc);
    rhoIn = kron(rhoh,kron(rhoc,rhow));
    dIn = diag(rhoIn);
    N.findSS(rhoIn);    
    rhoOut = N.rho;
    dOut = diag(rhoOut);
    s(j,:) = [sum(dIn([1:2:9 10:2:18 19:2:27])), sum(dOut([1:2:9 10:2:18 19:2:27])) sum(dIn([2:2:9 11:2:18 20:2:27])) sum(dOut([2:2:9 11:2:18 20:2:27]))];
    QC(1,j) = N.getTotalHeatFlow(2);
end
figure;hold on;
scatter(1:Nt,(QC(1,:)-qc)./qc);
title('Work,cold');

% r1 = [1 0; 0 1];
% N = model3SpinsFridgeTwoPhotonWC(g,T,k*E,E,'f','r','l',1);
% for j = 1:2
%     rhow = [r1(1,j)*(nw+1)/(2*nw+1) 1-r1(1,j) r1(1,j)*(nw)/(2*nw+1)];
%     rhoc = [r1(2,j)*(nc+1)/(2*nc+1) 1-r1(2,j) r1(2,j)*(nc)/(2*nc+1)];
%     rhow = diag(rhow);
%     rhoc = diag(rhoc);
%     N.findSS(kron(rhoh,kron(rhoc,rhow)));    
%     tmp = N.getTotalHeatFlow(2);
%     scatter(Nt+1,(tmp-qc)./qc,'r');
% end
% 
% % Two photons (HC)
% nh = 1./(exp(2*E(1)./T(1))-1);
% nc = 1./(exp(2*E(2)./T(2))-1);
% nw = 1./(exp(E(3)./T(3))-1);
% 
% rhow = exp(-E(3)./T(3)*(0:2));
% rhow = diag(rhow/sum(rhow));
% 
% Nt = 10000;
% r2 = rand(2,Nt);
% N = model3SpinsFridgeTwoPhotonHC(g,T,k*E,E,'f','r','l',1);
% for j = 1:Nt
%     rhoh = [r2(1,j)*(nh+1)/(2*nh+1) 1-r2(1,j) r2(1,j)*(nh)/(2*nh+1)];
%     rhoc = [r2(2,j)*(nc+1)/(2*nc+1) 1-r2(2,j) r2(2,j)*(nc)/(2*nc+1)];
%     rhoh = diag(rhoh);
%     rhoc = diag(rhoc);
%     N.findSS(kron(rhoh,kron(rhoc,rhow)));    
%     QC(2,j) = N.getTotalHeatFlow(2);
% end
% figure;
% scatter(1:Nt,(QC(2,:)-qc)./qc);
% title('Hot,cold');

% % Two photons (HW)
% nh = 1./(exp(2*E(1)./T(1))-1);
% nc = 1./(exp(E(2)./T(2))-1);
% nw = 1./(exp(2*E(3)./T(3))-1);
% 
% rhoc = exp(-E(2)./T(2)*(0:2));
% rhoc = diag(rhoc/sum(rhoc));
% 
% Nt = 10000;
% r3 = rand(2,Nt);
% N = model3SpinsFridgeTwoPhotonHW(g,T,k*E,E,'f','r','l',1);
% for j = 1:Nt
%     rhoh = [r3(1,j)*(nh+1)/(2*nh+1) 1-r3(1,j) r3(1,j)*(nh)/(2*nh+1)];
%     rhow = [r3(2,j)*(nw+1)/(2*nw+1) 1-r3(2,j) r3(2,j)*(nw)/(2*nw+1)];
%     rhoh = diag(rhoh);
%     rhow = diag(rhow);
%     N.findSS(kron(rhoh,kron(rhoc,rhow)));    
%     QC(3,j) = N.getTotalHeatFlow(2);
% end
% figure;
% scatter(1:Nt,(QC(3,:)-qc)./qc);
% title('Hot,work');
% 
% Nt = 1000;
% r4 = rand(3,Nt);
% N = model3SpinsFridgeTwoPhoton(g,T,k*E,E,'f','r','l',1);
% for j = 1:Nt
%     rhoh = [r4(1,j)*(nh+1)/(2*nh+1) 1-r4(1,j) r4(1,j)*(nh)/(2*nh+1)];
%     rhow = [r4(2,j)*(nw+1)/(2*nw+1) 1-r4(2,j) r4(2,j)*(nw)/(2*nw+1)];
%     rhoc = [r4(3,j)*(nc+1)/(2*nc+1) 1-r4(3,j) r4(3,j)*(nc)/(2*nc+1)];
%     rhoh = diag(rhoh);
%     rhow = diag(rhow);
%     rhoc = diag(rhoc);
%     N.findSS(kron(rhoh,kron(rhoc,rhow)));    
%     QC(4,j) = N.getTotalHeatFlow(2);
% end
% figure;
% scatter(1:Nt,(QC(4,:)-qc)./qc);
% title('Hot,work,cold');