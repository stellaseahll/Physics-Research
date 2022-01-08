load('optAngleF1.mat');
Fth = 1./(1+2*nbar).^2./(1+nbar)./nbar;
Fth = ones(500,1)*Fth;
% figure;

set(0,'defaulttextinterpreter','latex');
paper.width = 8;
paper.height = 8;
paper.FontSize = 20;
paper.LabelSize = 20;
figD = figure();
% paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
% paper.axheight = paper.axwidth;
% paper.height = paper.bottomgap+paper.axheight+paper.axgap;

% paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
% paper.xlabelOff = 0.09; %manual offset of xlabel from axis

% a = 0.9;
% fill([0 0 XEnd XEnd 0], [-0.1 0.8 0.8 -0.1 -0.1],[a a a],'linestyle','none');
theta(1,1) = 0;
theta(theta==max(max(theta))) = 0.055*pi;
division = 0:0.005:0.055;
contourf(nbar,gammat,theta/pi,division);
hold on;
plot(nbar,(1/3)./nbar,'r:','linewidth',3);
axis([0 10 0 1]);
set(gca,'YTick',[min(gammat) 0.5:0.5:1]);
set(gca,'XTick',[min(nbar) 5:5:10]);
set(gca,'XTickLabel',{'0','5','10'});
set(gca,'YTickLabel',{'0','0.5','1'});
xlabel('$\bar{n}$','interpreter','latex','FontSize',paper.LabelSize);
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',paper.LabelSize);
text(-2.1,0.95,'(a)','interpreter','latex','FontSize',20);
set(gca,'FontSize',paper.FontSize);
set(gca,'FontWeight','bold');
map = pink(length(division)-1);
map = map(end:-1:1,:);
colormap(map);
h=colorbar('eastoutside');
h.Label.String = '$\theta/\pi$';
h.Label.Interpreter = 'latex';
h.Label.FontSize = paper.LabelSize;
h.Ticks = 0:0.02:0.055;
h.FontSize = paper.FontSize;
% text(0.05,1.8,'(a)','interpreter','latex','FontSize',paper.LabelSize);
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca, 'Position', [pos(1) pos(2)+0.2 pos(3) pos(4)-0.15]);
% set(gca,'Position',pos);% set(figD.Children,'TickLabelInterpreter','latex');
set(figD.Children,'FontSize',paper.FontSize);
%print('fig2.eps','-dpsc2');
print('figF1angle.pdf','-dpdf');
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);



% a = 0.9;
% fill([0 0 XEnd XEnd 0], [-0.1 0.8 0.8 -0.1 -0.1],[a a a],'linestyle','none');
figD = figure();
F1 = -FI./F0;
division = 1:0.0005:1.0045;
F1(F1==max(max(F1))) = 1.0045;
contourf(nbar,gammat,F1,division); hold on;
plot(nbar,(1/3)./nbar,'r:','linewidth',2);
axis([0 10 0 1]);
set(gca,'YTick',[min(gammat) 0.5:0.5:1]);
set(gca,'XTick',[min(nbar) 5:5:10]);
set(gca,'XTickLabel',{'0','5','10'});
set(gca,'YTickLabel',{'0','0.5','1'});
xlabel('$\bar{n}$','interpreter','latex','FontSize',paper.LabelSize);
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',paper.LabelSize);
set(gca,'FontSize',paper.FontSize);
set(gca,'FontWeight','bold');
map = pink(length(division)-1);
map = map(end:-1:1,:);
colormap(map);
h=colorbar('eastoutside');
h.Label.String = '$\mathcal{F}_1/\mathcal{F}_1^{|g\rangle}$';
h.Label.Interpreter = 'latex';
h.Label.FontSize = paper.LabelSize;
h.Ticks = division(1:2:end);
h.FontSize = paper.FontSize;
% text(0.05,1.8,'(b)','interpreter','latex','FontSize',paper.LabelSize);
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca, 'Position', [pos(1) pos(2)+0.2 pos(3) pos(4)-0.15]);
% set(gca,'Position',pos);% set(figD.Children,'TickLabelInterpreter','latex');
set(figD.Children,'FontSize',paper.FontSize);
%print('fig2.eps','-dpsc2');
print('figF1_Fg.pdf','-dpdf');
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);
%
%
%



F1 = log10(FI./Fth);
F1(F1==min(min(F1)))= -3;
F1(F1==max(max(F1)))= 2;
division = -3:0.25:2;
[a b] = contour(nbar,gammat,F1,[-10 0 10]);
figD = figure();
contourf(nbar,gammat,F1,division); hold on;
plot(a(1,2:end),a(2,2:end),'b--','linewidth',3);
plot(nbar,(1/3)./nbar,'r:','linewidth',3);

axis([0 10 0 1]);
set(gca,'YTick',[min(gammat) 0.5:0.5:1]);
set(gca,'XTick',[min(nbar) 5:5:10]);
set(gca,'XTickLabel',{'0','5','10'});
set(gca,'YTickLabel',{'0','0.5','1'});
xlabel('$\bar{n}$','interpreter','latex','FontSize',paper.LabelSize);
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',paper.LabelSize);
text(-2.1,0.95,'(b)','interpreter','latex','FontSize',20);

set(gca,'FontSize',paper.FontSize);
set(gca,'FontWeight','bold');
map = pink(length(division)-1);
map = map(end:-1:1,:);
colormap(map);
h=colorbar('eastoutside');
h.Label.String = '$\log_{10} \mathcal{F}_1^{\rm opt}/\mathcal{F}_{\rm th}$';
h.Label.Interpreter = 'latex';
h.Label.FontSize = paper.LabelSize;
h.Ticks = division(1:4:end);
h.FontSize = paper.FontSize;
% text(0.05,1.8,'(b)','interpreter','latex','FontSize',paper.LabelSize);
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca, 'Position', [pos(1) pos(2)+0.2 pos(3) pos(4)-0.15]);
% set(gca,'Position',pos);% set(figD.Children,'TickLabelInterpreter','latex');
set(figD.Children,'FontSize',paper.FontSize);
%print('fig2.eps','-dpsc2');
print('figF1_Fth.pdf','-dpdf');
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);
