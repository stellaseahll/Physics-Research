set(0,'defaulttextinterpreter','latex');
paper.width = 8;
paper.height = 2.5;
paper.fontSize = 12;
c1 = [52 79 168]/255;
c2 = [193 24 0]/255;
c3 = [84 130 53]/255;
c4 = [207 138 68]/255;

x = load('compareModels.dat');
figD = figure();
semilogx(x(1,:),x(5,:)/0.01,'color', c1,'linewidth',2);
hold on;

semilogx(x(1,:),x(7,:)/0.01,'--','color',  c2,'linewidth',2);
semilogx(x(1,:),x(6,:)/0.01,':', 'color', c3,'linewidth',2);
fill([0.01 1 1 0.01 0.01],[0 0 -0.11 -0.11 0],[0.9 0.9 0.9],'EdgeColor','none'); hold on;
semilogx(x(1,:),x(2,:)/0.01,'color', c1,'linewidth',2); hold on;
hold on;

semilogx(x(1,:),x(4,:)/0.01,'--','color',  c2,'linewidth',2);
semilogx(x(1,:),x(3,:)/0.01,':', 'color', c3,'linewidth',2);

% 
text(0.1,-0.05,'\bf heating','FontSize',12,'FontWeight','bold','interpreter','latex');
% text(0.012,-0.08,'Heating','Rotation',90,'FontSize',12,'FontWeight','bold');
axis([0.01 1 -0.11 0.19]);
% grid on;
% set(gca,'FontSize', 18);
set(gca,'YTick',-0.1:0.1:0.1);
set(gca,'XTick',[0.01 0.1 1]);
set(gca,'XTickLabel',{0.01 0.1 1});
set(gca,'FontWeight','bold');
% set(gca,'XTick',0:1:5);
xlabel('$g/\Omega_{\mathrm {c}}$','interpreter','latex','FontSize',12);
ylabel('$\mathcal{J}_{\mathrm{c}}^\infty/\hbar\Omega_{\mathrm c}\kappa_{\mathrm c}$','interpreter','latex','FontSize',12);
set(gca,'FontWeight','bold');
% set(gca, 'XTickLabel', []);
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
set(figD.Children,'FontSize',paper.fontSize);
print('compareLGP2', '-dpdf', '-r600');
close(figD);
% export_fig heatHot.pdf;
