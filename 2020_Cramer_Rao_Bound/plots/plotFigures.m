set(0,'defaulttextinterpreter','latex');
paper.width = 8;
paper.height = 8;
paper.fontSize = 12;
load('ground_nbar1.541.mat');
% nbar = 1;
% filename = sprintf('XYZint_plus_nbar%.2f.mat',nbar);
% load(filename);
figD = figure();
plotContour(gt/pi,gammat,real(F5'/5/Fth));
set(gca,'YTick',0:2);
set(gca,'XTick',0:0.5:1);
xlabel('$g\tau_{\rm SA}/\pi$','interpreter','latex','FontSize',12);
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',12);
h=colorbar('eastoutside');
% mitte = (figD.Position(2)+figD.Position(4) + figD.Position(2))/2;
h.Label.String = '$F_5/5F_{\rm th}$';
h.Label.Interpreter = 'latex';
set(gca,'FontWeight','bold');
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);
filename = sprintf('groundState_F5_nbar%.2f.pdf',nbar);
print(filename, '-dpdf','-r600');
close(figD);
figD = figure();
plotContour(gt(2:99)/pi,gammat,X);
set(gca,'YTick',0:2);
set(gca,'XTick',0.1:0.4:0.9);
axis([0.1 0.9 0 2]);
xlabel('$g\tau_{\rm SA}/\pi$','interpreter','latex','FontSize',12);
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',12);
h=colorbar('eastoutside');
h.Ticks = 1:4;
h.Label.String = '$F_5/5F_1$';
h.Label.Interpreter = 'latex';
set(gca,'FontWeight','bold');
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);
filename = sprintf('groundState_F5F1_nbar%.2f.pdf',nbar);
print(filename, '-dpdf','-r600');
close(figD);