% Time evolutions for resonant exchange thermalization
% Paper plots

%% Testing part

J = 1/2;
omegaP = 10;
g = 1; 
%assuming dipole-dipole coupling (or purely resonant exchange)
gsb = diag(g*[1,0,0,0,0]);
gamma = 1e-2;

T = 15;

rhoth = exp(-omegaP/T*(0:2*J)).';
rhoth = rhoth/sum(rhoth);

rho0 = zeros(2*J+1,1);
rho0(1)=1;

%varying parameters
omegaS = omegaP;
tau1 = 2*pi/(g/2);
Omega = omegaS+omegaP;
tau2 = 2*pi/(sqrt((g/2)^2+Omega^2)/2);
tau = tau1*(1);
% ps = {'-'};
% pc = {'b','k','m','r','y'};

Nt = 2001;
t = logspace(0,10,Nt);

%work powers VS t
W = cell( length(omegaS),length(tau) );
%trace distances to P-thermal VS t
D = cell( length(omegaS),length(tau) );

for n1=1:length(omegaS)
    rhothS = exp(-omegaP/T*(0:2*J)).';
    rhothS = rhothS/sum(rhothS);
    
    for n2=1:length(tau)
        tic
        M = collisionModelSpinBath(2*J+1,omegaS(n1),2,omegaP,T,gamma,tau(n2),gsb,diag(rho0));
        M.prepareSim();
        [rho, obst, w] = M.runMEScatter(t,{diag(rhoth)}); %compute overlap w/ thermal state also
        obs{n1,n2} = obst;
        W{n1,n2} = cumtrapz(t,w);
        %trace distance to thermal state
        D{n1,n2} = squeeze( sum(sum(rho.*rho,1),2) ).' - 2*obst + sum(rhoth.^2);
        toc
    end
    
end

%Test Plot

% figure();
% semilogx(t,real(D{1,1}),ps{1});
% hold on
% 
% for n1=1:length(omegaS)
%     for n2=1:length(tau)
%         if ~(n1==1 && n2==1)
%             semilogx(t,real(D{n1,n2}),[pc{n1},ps{n2}]);
%         end
%     end
% end
% 
% hold off
% 
figure();
% semilogx(t,real(W{1,1}),ps{1});
hold on
for n1=1:length(omegaS)
    for n2=1:length(tau)
%         if ~(n1==1 && n2==1)
            semilogx(t,real(W{n1,n2}));
%         end
    end
end

hold off
% 
% %% Produce data for paper plots!
% J = 2;
% omegaP = 10;
% g = 1; 
% %assuming dipole-dipole coupling (or purely resonant exchange)
% gsb = diag(g*[1,1,1,0,0]);
% gamma = 1e-2;
% T = 15;
% Nt = 100001;
% 
% rhoth = exp(-omegaP/T*(0:2*J)).';
% rhoth = rhoth/sum(rhoth);
% 
% rho0 = zeros(2*J+1,1);
% rho0(1)=1;
% 
% 
% %varying parameters
% omegaS = omegaP*[0.8,0.9,1.1,1.2];
% tau = [1e-1,10];
% dt = [5, 5, 5, 5; 1, 1, 1, 1].';
% dt = [5, 5, 5, 5; 5, 5, 5, 5].';
% ps = {'--','-'};
% pc = {'g','b','k','r'};
% 
% %work powers and cum work VS t
% PW = cell( length(omegaS),length(tau) );
% W = cell( length(omegaS),length(tau) );
% PW2 = cell( length(omegaS),length(tau) );
% W2 = cell( length(omegaS),length(tau) );
% ts = cell( length(omegaS),length(tau) );
% 
% for n1=1:length(omegaS)
%     for n2=1:length(tau)
%         tic
%         t = (0:Nt)*dt(n1,n2);
%         M = collisionModelSpinBath(2*J+1,omegaS(n1),2,omegaP,T,gamma,tau(n2),gsb,diag(rho0));
%         M.prepareSim();
%         [rho, obst, w] = M.runMEScatter(t,{diag(rhoth)}); %compute overlap w/ thermal state also
%         PW{n1,n2} = w;
%         W{n1,n2} = cumtrapz(t,w);
%         M = collisionModelSpinBath(2*J+1,omegaS(n1),2,omegaP,T,gamma,tau(n2)*3,gsb/3,diag(rho0));
%         M.prepareSim();
%         [rho, obst, w] = M.runMEScatter(t,{diag(rhoth)}); %compute overlap w/ thermal state also
%         PW2{n1,n2} = w;
%         W2{n1,n2} = cumtrapz(t,w);
%         ts{n1,n2} = t;
%         toc
%     end
%     
% end
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
% %% Paper Plots a, b
% 
% paper.width = 8;
% paper.height = 4;
% paper.fontSize = 8;
% set(0,'defaulttextinterpreter','latex')
% 
% c = colormap(lines);
% 
% figD = figure(); 
% 
% semilogx(gamma*[10, 1e6],[0, 0],'Color',0.33*[1 1 1],'LineWidth',0.5);
% hold on
% semilogx(gamma*ts{1,1},real(W{1,1})/omegaP,'Color',c(1,:),'LineWidth',1);
% semilogx(gamma*ts{1,2},real(W{1,2})/omegaP,':','Color',c(1,:),'LineWidth',1.5);
% 
% semilogx(gamma*ts{2,1},real(W{2,1})/omegaP,'Color',c(4,:),'LineWidth',1);
% semilogx(gamma*ts{2,2},real(W{2,2})/omegaP,':','Color',c(4,:),'LineWidth',1.5);
% 
% semilogx(gamma*ts{3,1},real(W{3,1})/omegaP,'Color',c(2,:),'LineWidth',1);
% semilogx(gamma*ts{3,2},real(W{3,2})/omegaP,':','Color',c(2,:),'LineWidth',1.5);
% 
% semilogx(gamma*ts{4,1},real(W{4,1})/omegaP,'Color',c(3,:),'LineWidth',1);
% semilogx(gamma*ts{4,2},real(W{4,2})/omegaP,':','Color',c(3,:),'LineWidth',1.5);
% 
% 
% hold off
% xlim(gamma*[20,4e5]);
% ylim([-2,2]/omegaP);
% 
% set(gca,'FontWeight','bold')
% set(gca,'XTick',[100,1e3,1e4,1e5]*gamma)
% set(gca,'XTickLabel',{1,10,100,1000})
% %set(gca,'YTick',[1,2,3,4,5,6])
% xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
% ylabel('\bf work $W/ \hbar \omega_p$','interpreter','latex','FontSize',8)
% 
% text(gamma*30,1.5/omegaP,'\bf (a)','interpreter','latex','FontSize',9)
% 
% %for saving
% figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'FontSize',paper.fontSize)
% figD.PaperPositionMode = 'auto';
% posD=get(figD,'Position');
% set(figD,'PaperSize',[posD(3), posD(4)])
% set(figD, 'color', 'none');
% %print('fig1', '-dpdf', '-r600')
% export_fig fig1a.pdf
% close(figD)
% 
% 
% figD = figure(); 
% 
% semilogx(gamma*[10, 1e6],[0, 0],'Color',0.33*[1 1 1],'LineWidth',0.5);
% hold on
% semilogx(gamma*ts{1,1},real(W2{1,1})/omegaP,'Color',c(1,:),'LineWidth',1);
% semilogx(gamma*ts{1,2},real(W2{1,2})/omegaP,':','Color',c(1,:),'LineWidth',1.5);
% 
% semilogx(gamma*ts{2,1},real(W2{2,1})/omegaP,'Color',c(4,:),'LineWidth',1);
% semilogx(gamma*ts{2,2},real(W2{2,2})/omegaP,':','Color',c(4,:),'LineWidth',1.5);
% 
% semilogx(gamma*ts{3,1},real(W2{3,1})/omegaP,'Color',c(2,:),'LineWidth',1);
% semilogx(gamma*ts{3,2},real(W2{3,2})/omegaP,':','Color',c(2,:),'LineWidth',1.5);
% 
% semilogx(gamma*ts{4,1},real(W2{4,1})/omegaP,'Color',c(3,:),'LineWidth',1);
% semilogx(gamma*ts{4,2},real(W2{4,2})/omegaP,':','Color',c(3,:),'LineWidth',1.5);
% 
% 
% hold off
% xlim(gamma*[20,4e5]);
% ylim([-2,2]/omegaP);
% 
% set(gca,'FontWeight','bold')
% set(gca,'XTick',[100,1e3,1e4,1e5]*gamma)
% set(gca,'XTickLabel',{1,10,100,1000})
% %set(gca,'YTick',[1,2,3,4,5,6])
% xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
% ylabel('\bf work $W/ \hbar \omega_p$','interpreter','latex','FontSize',8)
% 
% text(gamma*30,1.5/omegaP,'\bf (b)','interpreter','latex','FontSize',9)
% 
% %for saving
% figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'FontSize',paper.fontSize)
% figD.PaperPositionMode = 'auto';
% posD=get(figD,'Position');
% set(figD,'PaperSize',[posD(3), posD(4)])
% set(figD, 'color', 'none');
% %print('fig1', '-dpdf', '-r600')
% export_fig fig1b.pdf
% close(figD)