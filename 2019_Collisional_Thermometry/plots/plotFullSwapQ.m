clc;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
% c = colormap(lines);
figD = figure(); 
c = [0.6:-0.1:0;0.6:-0.1:0;0.6:-0.1:0]
% close(figure);

%% Plot Fisher Info for different K for 1 pass
% figD = figure();  
G = 0.5;
gt = 1.0;
nbar = 1;
np = 7;
isSteady = 1;
filename = sprintf('G%.3fpi,gt%.3f,nbar%.3f,np%d,isSteady%d.mat',G,gt,nbar,np,isSteady);
load(filename);
%% 
for i = 1:np
    plot(theta,real(F(:,i)),'color',c(:,i),'LineStyle','-','LineWidth',1); hold on;
end
% plot(theta,real(Fz(:,1)),'color',c(1,:),'LineStyle','--','LineWidth',1);
% plot(theta,real(F(:,2)),'color',c(2,:),'LineStyle','-','LineWidth',1);
% plot(theta,real(Fz(:,2)),'color',c(2,:),'LineStyle','--','LineWidth',1);
% plot(theta,real(F(:,3)),'color',c(3,:),'LineStyle','-','LineWidth',1);
% plot(theta,real(Fz(:,3)),'color',c(3,:),'LineStyle','--','LineWidth',1);

set(gca,'FontWeight','bold');
xlabel('\bf initial state $\phi/\pi $','interpreter','latex','FontSize',8);
ylabel('\bf $F_N$','interpreter','latex','FontSize',8);
% xtxt = 1.8;
% ytxt = 0.86;
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
print('fig4', '-dpdf', '-r600')
close(figD);
