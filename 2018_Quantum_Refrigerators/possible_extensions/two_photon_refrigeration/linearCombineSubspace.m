k = 0.001;
g = logspace(-3,0,100);
g = g(1:end-1);
E = [5 1 4];
TW = 5;
TH = 2;
TC = 1;
T = [TH TC TW];
r = 0:0.025:1;
nh = 1./(exp(E(1)./T(1))-1);
nc = 1./(exp(2*E(2)./T(2))-1);
nw = 1./(exp(2*E(3)./T(3))-1);
rhoh = exp(-E(1)./T(1)*(0:2));
rhoh = diag(rhoh/sum(rhoh));
for i = 1:length(g)
%     % One photon
%     M = model3SpinsFridge(g,T,k*E,E,'f','r','l',1);
%     M.findSS();
%     qc = M.getTotalHeatFlow(2);
    i
    % Two photons (WC)
    Nl = model3SpinsFridgeTwoPhotonWC(g(i),T,k*E,E,'o','r','l',1);
    Ng = model3SpinsFridgeTwoPhotonWC(g(i),T,k*E,E,'o','r','g',1);
    Np = model3SpinsFridgeTwoPhotonWC(g(i),T,k*E,E,'o','r','p',1);
    rhow = diag([1 0 0]);
    rhoc = diag([1 0 0]);
    rho1 = kron(rhoc,rhow);

    rhow = diag([1 0 0]);
    rhoc = diag([0 1 0]);
    rho2 = kron(rhoc,rhow);
    for j = 1:length(r)
        rhocw = r(j)*rho1 + (1-r(j))*rho2;
        rhocw = rhocw/trace(rhocw);
        Nl.findSS(kron(rhoh,rhocw));   
        Ng.findSS(kron(rhoh,rhocw));   
        Np.findSS(kron(rhoh,rhocw));   
        QCl(i,j) = Nl.getTotalHeatFlow(2);
        QCg(i,j) = Ng.getTotalHeatFlow(2);
        QCp(i,j) = Np.getTotalHeatFlow(2);
    end
    save 2photonWC_lgp_Varyg_TC1_TH2_TW5_EH5_EW4_EC1.mat;
%     QCrat(i,:) = (QC(i,:)-qc)/qc;
end
% save linearCombineSubspace_VaryTW_TC1_TH2_EH5_EW4_EC1.mat;
% load linearCombineSubspace_VaryTW_TC1_TH2_EH5_EW4_EC1.mat;
% hold on;
% figure;
% plot(r,QCrat);
% legend('proportion of rho1','QC-qc/qc')
% figure;
% plot(TW,QCrat(:,1))
% legend('TW','QC-qc/qc')

% k = [0 1];
% k = k/max(k);
% map = [red;1 1 1];
% QCrat(end,end) = -0.8;
% map = interp1(k,map,linspace(0,1,length(-0.8:0.1:0)));
% map = [map(1:end-1,:);blue];
% contourf(r,TW,QCrat,-0.8:0.1:0.1)
% colormap(map);
% ylabel('Work Bath Temperature');
% set(gca,'FontSize', 18);
% print('test.eps','-dpsc2');

figure;
semilogx(g,QCg(:,end),'b:','linewidth',2);
hold on;
semilogx(g,QCl(:,end),'b--','linewidth',2);
semilogx(g,QCp(:,end),'b','linewidth',2);
semilogx(g,QCg(:,1),'r:','linewidth',2);
semilogx(g,QCl(:,1),'r--','linewidth',2);
semilogx(g,QCp(:,1),'r','linewidth',2);
