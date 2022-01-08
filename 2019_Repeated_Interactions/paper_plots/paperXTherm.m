% Produce data for paper plots
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

rho0 = zeros(2*J+1,1);
rho0(1)=1;

% tau = [1e-1,10];
% dt = [5, 5, 5, 5; 1, 1, 1, 1].';
% dt = [5, 5, 5, 5; 5, 5, 5, 5].';
% ps = {'--','-'};
% pc = {'g','b','k','r'};

%work powers and cum work VS t
dtau = 2*pi/(g/2); % period of oscillation
tau = [logspace(log10(dtau/2000),log10(dtau/20),300),dtau/20:dtau/100:2.5*dtau dtau*([0.9:0.00002:1.1 1.9:0.00002:2.1])];
G = g/2;
Omega = omegaS+omegaP;
dtau = 2*pi/(sqrt(G^2+Omega^2)/2);
tau = [tau dtau/20:dtau/100:2.5*dtau/g*G];
tau = unique(tau);
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
        p0(:,n) = diag(rho);
        Ez(n) = sum(p0(:,n).*Jz);
        Erg(n) = Ez(n) - sum(sort(p0(:,n),'descend').*Jz);
end
Erg = real(Erg)>0;
Erg = Erg*1.0;
Erg(Erg==0) = NaN;

%% Paper Plots b

paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')
% 
c = colormap(lines);
xend = 30;
y0 = min([real(Ez) Ezth])*1.05;
y1 = max(real(Ez))*1.1;

figD = figure(); 
idx = find(Erg>0);
k = diff(idx);
idx2 = [idx(1) idx(find(k~=1)+1); idx(find(k~=1)) idx(end)];
hold on;
for i = 1:length(idx2)
    
    fill([g*tau(idx2(1,i)) g*tau(idx2(1,i)) g*tau(idx2(2,i)) g*tau(idx2(2,i)) g*tau(idx2(1,i))],[y0 y1 y1 y0 y0],c(3,:),'EdgeColor','none');
end
plot(g*tau,Ez,'color',c(1,:),'LineWidth',1.0);
% plot(g*tau(Erg>0),Ez(Erg>0),'Color',c(2,:),'LineWidth',1.0);
% plot([g*tau g*tau],[Erg*(y1-(y1-y0)*0.01) Erg*(y0)],'Color',[1 0 0],'LineWidth',2);
plot([0 g*tau(end)],[Ezth Ezth],'Color',c(2,:),'LineStyle','--','LineWidth',0.5);
% plot([0 g*tau(end)],[E1 E1],'Color',0.33*[1 1 1],'LineStyle','--','LineWidth',0.33);

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
hold off
xlim([0,30]);

ylim([y0 y1]);
% 
set(gca,'FontWeight','bold')
set(gca,'XTick',0:5:xend)
set(gca,'XTickLabel',{0,5,10,15,20,25,30})
box on
%set(gca,'YTick',[1,2,3,4,5,6])
xlabel('\bf interaction time $g_x \tau$','interpreter','latex','FontSize',8)
ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% 
text(0.5,y1 - 0.1*(y1-y0),'\bf (a)','interpreter','latex','FontSize',9)
% 
% for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig2a_g25', '-dpdf', '-r600')
% export_fig fig2a.pdf
close(figD)

J = 2;


rhoth = exp(-omegaP/T*(0:2*J)).';
rhoth = rhoth/sum(rhoth);

rho0 = zeros(2*J+1,1);
rho0(1)=1;

% tau = [1e-1,10];
% dt = [5, 5, 5, 5; 1, 1, 1, 1].';
% dt = [5, 5, 5, 5; 5, 5, 5, 5].';
% ps = {'--','-'};
% pc = {'g','b','k','r'};

%work powers and cum work VS t
dtau = 2*pi/(g/2); % period of oscillation
tau = [logspace(log10(dtau/2000),log10(dtau/20),300),dtau/20:dtau/100:2.5*dtau dtau*([0.9:0.00002:1.1 1.9:0.00002:2.1])];
G = g/2;
Omega = omegaS+omegaP;
dtau = 2*pi/(sqrt(G^2+Omega^2)/2);
tau = [tau dtau/20:dtau/100:2.5*dtau/g*G];
tau = unique(tau);
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
        p0(:,n) = diag(rho);
        Ez(n) = sum(p0(:,n).*Jz);
        Erg(n) = Ez(n) - sum(sort(p0(:,n),'descend').*Jz);
end
Erg = real(Erg)>0;
Erg = Erg*1.0;
Erg(Erg==0) = NaN;

%% Paper Plots b

paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')
% 
c = colormap(lines);
xend = 30;
y0 = min([real(Ez) Ezth])*1.05;
y1 = max(real(Ez))*1.1;

figD = figure(); 
idx = find(Erg>0);
k = diff(idx);
idx2 = [idx(1) idx(find(k~=1)+1); idx(find(k~=1)) idx(end)];
hold on;
for i = 1:length(idx2)
    
    fill([g*tau(idx2(1,i)) g*tau(idx2(1,i)) g*tau(idx2(2,i)) g*tau(idx2(2,i)) g*tau(idx2(1,i))],[y0 y1 y1 y0 y0],c(3,:),'EdgeColor','none');
end
plot(g*tau,Ez,'color',c(1,:),'LineWidth',1.0);
% plot(g*tau(Erg>0),Ez(Erg>0),'Color',c(2,:),'LineWidth',1.0);
% plot([g*tau g*tau],[Erg*(y1-(y1-y0)*0.01) Erg*(y0)],'Color',[1 0 0],'LineWidth',2);
plot([0 g*tau(end)],[Ezth Ezth],'Color',c(2,:),'LineStyle','--','LineWidth',0.5);
% plot([0 g*tau(end)],[E1 E1],'Color',0.33*[1 1 1],'LineStyle','--','LineWidth',0.33);

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
hold off
xlim([0,30]);

ylim([y0 y1]);
% 
set(gca,'FontWeight','bold')
set(gca,'XTick',0:5:xend)
set(gca,'XTickLabel',{0,5,10,15,20,25,30})
box on
%set(gca,'YTick',[1,2,3,4,5,6])
xlabel('\bf interaction time $g_x \tau$','interpreter','latex','FontSize',8)
ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% 
text(0.5,y1 - 0.1*(y1-y0),'\bf (b)','interpreter','latex','FontSize',9)
% 
% for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig2b_g25', '-dpdf', '-r600')
% export_fig fig2a.pdf
close(figD)