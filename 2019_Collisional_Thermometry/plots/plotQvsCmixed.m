clc;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure(); 

% close(figure);

%% Plot Fisher Info for different K for 1 pass
% figD = figure();  
G = 0.01;
gt = 1.0;
nbar = 1;
np = 7;
isSteady = 1;
filename = sprintf('DiffMixedState,G%.3fpi,gt%.3f,nbar%.3f,np%d,isSteady%d.mat',G,gt,nbar,np,isSteady);
load(filename);
%% 
plot(z,real(F(:,2)),'color',c(2,:),'LineStyle','-','LineWidth',1);hold on;
plot(z,real(Fz(:,2)),'color',c(2,:),'LineStyle','--','LineWidth',1);
plot(z,real(F(:,1)),'color',c(1,:),'LineStyle','-','LineWidth',1); hold on;

plot(z,real(Fz(:,1)),'color',c(1,:),'LineStyle','--','LineWidth',1);

% plot(theta,real(F(:,3)),'color',c(3,:),'LineStyle','-','LineWidth',1);
% plot(theta,real(Fz(:,3)),'color',c(3,:),'LineStyle','--','LineWidth',1);

set(gca,'FontWeight','bold');
xlabel('initial state \bf $p_A $','interpreter','latex','FontSize',8);
ylabel('\bf $F_N$','interpreter','latex','FontSize',8);
xtxt = 1.8;
ytxt = 0.86;
% text(xtxt,max(F(:,2))*ytxt,'\bf (a)','interpreter','latex','FontSize',9)


%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('figMixed1a', '-dpdf', '-r600')
close(figD);
%%
figD = figure(); 
r = alpha(:,1);
plot(z,real(r),'color',c(1,:),'LineStyle','-','LineWidth',1); hold on;
plot(z,real(alphaz(:,1)),'color',c(1,:),'LineStyle','--','LineWidth',1);
r = alpha(:,end);
plot(z,real(r),'color',c(2,:),'LineStyle','-','LineWidth',1); hold on;
plot(z,real(alphaz(:,end)),'color',c(1,:),'LineStyle','--','LineWidth',1);
set(gca,'FontWeight','bold');
xlabel('initial state \bf $p_A $','interpreter','latex','FontSize',8);
ylabel('\bf $\alpha_N$','interpreter','latex','FontSize',8);
% text(xtxt,max(r)*ytxt,'\bf (b)','interpreter','latex','FontSize',9)
axis([0 1 0.9 1.1])
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
% print('figMixed3b', '-dpdf', '-r600')
close(figD);