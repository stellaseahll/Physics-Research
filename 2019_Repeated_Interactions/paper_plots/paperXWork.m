%% Produce data for paper plots
% clear;clc;clf;
figure;
J = 0.5;
omegaP = 11;
omegaS = 11;
g = 2.5; 
%assuming dipole-dipole coupling (or purely resonant exchange)
gsb = diag(g*[1,0,0,0,0]);
gamma = g/1000;
T = 1.5*omegaP;
Nt = 100001;

rhoth = exp(-omegaP/T*(0:2*J)).';
rhoth = rhoth/sum(rhoth);
eta = exp(-omegaP/T*(0:1)).';
eta = eta/sum(eta);
eta = diag(eta);
rho0 = zeros(2*J+1,1);
rho0(1)=1;

% tau = [1e-1,10];
% dt = [5, 5, 5, 5; 1, 1, 1, 1].';
% dt = [5, 5, 5, 5; 5, 5, 5, 5].';
% ps = {'--','-'};
% pc = {'g','b','k','r'};

%work powers and cum work VS t
dtau = 2*pi/(g/2); % period of oscillation
G = g/2;
Omega = omegaS+omegaP;
% dtau = 2*pi/(sqrt(G^2+Omega^2)/2);
tau = [0.75:0.0005:1.25]*dtau;%[dtau/100:dtau/100:10*dtau];%
% tau = [logspace(log10(dtau/2000),log10(dtau/20),100),dtau/20:dtau/2000:2.2*dtau];
% tau = [0.1:dtau/400:dtau*2];
Jz = (-J:J)';
Ezth = sum(Jz.*rhoth);
p0 = zeros(2*J+1,length(tau));
Ez = zeros(1,length(tau));
Erg = zeros(1,length(tau));
% PW = cell( length(omegaS),length(tau) );
% W = cell( length(omegaS),length(tau) );
% PW2 = cell( length(omegaS),length(tau) );
% W2 = cell( length(omegaS),length(tau) );
% ts = cell( length(omegaS),length(tau) );

for n=1:length(tau)
        M = collisionModelSpinBath(2*J+1,omegaS,2,omegaP,T,gamma,tau(n),gsb,diag(rho0));
        M.findSSwithoutSim();
        rho = M.rhoSSscatter{1};
        Rho = kron(eta,rho);
        UL = expm(-1i*M.H0{1}*tau(n)/2);
        UR = UL';
        E = M.S{1}'*M.H0{1}*M.S{1}-M.H0{1};
        W(n) = sum(sum((UL*Rho*UR).*E));
        U = expm(-1i*(M.Hint{1} + M.H0{1})*tau(n));
        Comm(n) = norm(U*M.H0{1} - M.H0{1}*U);
        p0(:,n) = diag(rho);
        Ez(n) = sum(p0(:,n).*Jz);
        Erg(n) = Ez(n) - sum(sort(p0(:,n),'descend').*Jz);
end
Erg = real(Erg)>0;
Erg = Erg*1.0;
Erg(Erg==0) = NaN;
% 
% save('workResonantExchangePlot.mat');
% 
% figure();
% semilogx(ts{1,1},real(W{1,1}),[pc{1},ps{1}]);
% grid on
% hold on
% for n1=1:length(omegaS)
%     for n2=1:length(tau)
%         if ~(n1==1 && n2==1)
%             semilogx(ts{n1,n2},real(W{n1,n2}),[pc{n1},ps{n2}]);
%         end
%     end
% end
% hold off
% 
% 
% 
%% Paper Plots b

paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

% 
c = colormap(lines);
xstart = 10;
xend = 15;
y0 = 0;
y1 = 0.0022;

figD = figure(); 
set(figD,'defaultAxesColorOrder',[c(2,:);c(1,:)]);
yyaxis left
% plot(g*tau,W/omegaP,'color',c(1,:),'linewidth',1);
idx = find(g*tau<=15);
plot(g*tau(idx),W(idx)/omegaP,'color',c(2,:),'linewidth',1);
ylim([y0 y1]);

set(gca,'YTick',0:0.001:0.003);
ylabel('\bf power $\dot{W}/\hbar\gamma\omega_p$','interpreter','latex','FontSize',8)
yyaxis right
plot(g*tau,Ez,'color',c(1,:),'linewidth',1,'linestyle','--'); hold on;
set(gca,'XTick',10:15)
set(gca,'YTick',-0.15:0.1:0.15);

ylim([Ezth -Ezth]);
xlim([xstart xend]);
ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)

% idx = find(W>0);
% k = diff(idx);
% idx2 = [idx(1) idx(find(k~=1)+1); idx(find(k~=1)) idx(end)];
% hold on;
% for i = 1:length(idx2)
%     
%     fill([g*tau(idx2(1,i)) g*tau(idx2(1,i)) g*tau(idx2(2,i)) g*tau(idx2(2,i)) g*tau(idx2(1,i))],[y0 y1 y1 y0 y0],c(1,:),'EdgeColor','none');
% end
% plot(g*tau,W,'k','LineWidth',1.0);
% plot([0 g*tau(end)],[Ezth Ezth],'Color',0.01*[1 1 1],'LineStyle','--','LineWidth',0.5);
% hold off
% xlim([0,g*tau(end)]);
% 
% ylim([y0 y1]);
% % 
set(gca,'FontWeight','bold')
% set(gca,'XTick',10:15)
% set(gca,'XTickLabel',{10,15,20,25})
box on
%set(gca,'YTick',[1,2,3,4,5,6])
xlabel('\bf interaction time $g_x\tau$','interpreter','latex','FontSize',8)
% ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% % 
% text(0.5,y1 - 0.1*(y1-y0),'\bf (b)','interpreter','latex','FontSize',9)
% % 
% % for saving
% figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'FontSize',paper.fontSize)
% figD.PaperPositionMode = 'auto';
% posD=get(figD,'Position');
% set(figD,'PaperSize',[posD(3), posD(4)])
% set(figD, 'color', 'none');
% print('fig3_color', '-dpdf', '-r600')
% % export_fig fig2a.pdf
% close(figD)