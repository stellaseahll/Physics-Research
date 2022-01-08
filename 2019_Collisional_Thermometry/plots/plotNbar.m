clc;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
G = 0.25;
gt = 0.1;
% nbar = 1;
np = 7;
isSteady = 1;
Theta = 0.5;
figD = figure();
for idxx = 1:length(Theta)
    filename = sprintf('DiffNbar,G%.3fpi,gt%.3f,theta%.3fpi,np%d,isSteady%d.mat',G,gt,Theta(idxx),np,isSteady);
    load(filename);
    semilogx(nbar,Fz(:,end),'color',c(idxx,:),'LineStyle',':','LineWidth',3);hold on;
    semilogx(nbar,F(:,end),'color',c(idxx,:),'LineStyle','-','LineWidth',2);hold on;
end
set(gca,'FontWeight','bold');
xlabel('\bf mean occ. num. $\bar{n} $','interpreter','latex','FontSize',8);
ylabel('\bf $F_7$','interpreter','latex','FontSize',8);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('figNbarPlus', '-dpdf', '-r600');
close(figD);
% 
% paper.width = 8;
% paper.height = 3;
% paper.fontSize = 8;
% figD = figure();
% figD = figure();
% for idxx = 1:length(Theta)
%     filename = sprintf('DiffNbar,G%.3fpi,gt%.3f,theta%.3fpi,np%d,isSteady%d.mat',G,gt,Theta(idxx),np,isSteady);
%     load(filename);
%     for j = 1:length(F)
%         alpha(j,:) = M.getAllExponent(F(j,:));
%         alphaz(j,:) = M.getAllExponent(Fz(j,:));
%     end
%     semilogx(nbar,alpha(:,end),'color',c(idxx,:),'LineStyle','-','LineWidth',2);hold on;
%     semilogx(nbar,alphaz(:,end),'color',c(idxx,:),'LineStyle',':','LineWidth',3);hold on; 
% end
% set(gca,'FontWeight','bold');
% xlabel('\bf mean occ. num. $\bar{n} $','interpreter','latex','FontSize',8);
% ylabel('\bf $\alpha_7$','interpreter','latex','FontSize',8);
% figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'FontSize',paper.fontSize);
% figD.PaperPositionMode = 'auto';
% posD=get(figD,'Position');
% set(figD,'PaperSize',[posD(3), posD(4)]);
% set(figD, 'color', 'none');
% print('figNbarExc', '-dpdf', '-r600');
% close(figD);