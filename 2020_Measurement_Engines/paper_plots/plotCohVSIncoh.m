%% load data
panels = {'(a)','(b)'};

%load('cohVSincoh_rates_wp1_ws100_g3.54_kh0.0010_kc0.1000_nc1.00_lowDim.mat'); 
load('cohVSincoh_rates_wp1_ws100_g3.54_kh0.0010_kc0.1000_nc1.00.mat'); 
wi = real(-WI)/kh/ws; ni = real(-WI)./real(QhI); 
wc = real(-WC)/kh/ws; nc = real(-WC)./real(QhC); 

%% plot wide fig (small g) 

paper.leftgap = 1.3;
paper.bottomgap = 0.9;
paper.axgap = 0.5;
paper.width = 14;
paper.height = 6;
paper.N = 2;
set(0,'defaulttextinterpreter','latex');
figD = figure(666);
c = colormap(lines);
% paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
% paper.axheight = paper.axwidth;
% paper.height = paper.bottomgap+paper.axheight+paper.axgap;

paper.fontSize = 12;
paper.labelSize = 12;
% paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
% paper.xlabelOff = 0.09; %manual offset of xlabel from axis

ax=subplot(1,2,1); 
% a = 0.9;
% fill([0 0 XEnd XEnd 0], [-0.1 0.8 0.8 -0.1 -0.1],[a a a],'linestyle','none');
semilogx(km/wp,wi,'Color','k','LineWidth',1);hold on;
semilogx(G/wp,wc,'Color',c(1,:),'LineWidth',1);
xlabel('$\gamma / \omega$, $\zeta / \omega$','interpreter','latex','FontSize',paper.labelSize);
ylabel('Power $\dot{W}/\hbar\Omega\kappa_h$','interpreter','latex','FontSize',paper.labelSize);
xlim([1e-4 1]);
ylim([0 0.8]);
%set(ax,'XTick',-8:2:8);
%set(ax,'YTick',0:0.1:0.4);
set(ax,'XTick',[1e-4,1e-3,1e-2,1e-1,1],'XTickLabels',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','1'});
% set(ax,'YTick',0:0.2:0.8);
text(2e-4,0.8*0.87,'(a)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);

ax = subplot(1,2,2);
semilogx(km/wp,ni,'Color','k','LineWidth',1);hold on;
semilogx(G/wp,nc,'Color',c(1,:),'LineWidth',1);
xlabel('$\gamma / \omega$, $\zeta / \omega$','interpreter','latex','FontSize',paper.labelSize);
ylabel('Efficiency $\eta$','interpreter','latex','FontSize',paper.labelSize);

xlim([1e-4 1]);
ylim([0 1]);
set(ax,'XTick',[1e-4,1e-3,1e-2,1e-1,1],'XTickLabels',{'$10^{-4}$','$10^{-3}$','$10^{-2}$','$10^{-1}$','1'});
%set(ax,'XTick',-8:2:8);
set(ax,'YTick',0:0.2:1);
text(2e-4, 0.87,'(b)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);

% 
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'TickLabelInterpreter','latex');
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
% %print('fig2.eps','-dpsc2');
print('figure4.pdf', '-dpdf', '-r600')
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);