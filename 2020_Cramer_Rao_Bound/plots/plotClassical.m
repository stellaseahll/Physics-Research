clear;
%% load data
panels = {'(a)','(b)'};

% load('8AncillaT2_g_e_p_tol6.mat');
% Fg = Fs{1}(1,:);
% Fe = Fs{2}(1,:);
% Fp = Fs{3}(1,:);

%here saved as arrays for a fixed gt=pi/100
load('ZZint_plus_nbar1.54.mat');
load('Classical_nbar1.541.mat');

X = real(FN./F1)'/5;
X(X<1) = 1;
%% plot wide fig (small g)

paper.leftgap = 1.3;
paper.bottomgap = 0.9;
paper.axgap = 0.5;
paper.width = 14;
paper.height = 3.8;
paper.N = 2;
set(0,'defaulttextinterpreter','latex');
figD = figure(666);

% paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
% paper.axheight = paper.axwidth;
% paper.height = paper.bottomgap+paper.axheight+paper.axgap;
paper.fontSize = 30;
paper.labelSize = 36;
% paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
% paper.xlabelOff = 0.09; %manual offset of xlabel from axis

ax=subplot(1,2,1);
% a = 0.9;
% fill([0 0 XEnd XEnd 0], [-0.1 0.8 0.8 -0.1 -0.1],[a a a],'linestyle','none');
map=plotContour(gt/pi,gammat,real(FN'/5/Fth));
set(gca,'YTick',0:2);
set(gca,'XTick',0:0.5:1);
xlabel('$g\tau_{\rm SA}/\pi$','interpreter','latex','FontSize',paper.labelSize );
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',paper.labelSize );
axis([0 1 0 2]);
colormap(gca,map);
h=colorbar('eastoutside');
h.Ticks = 0:1:4;
% mitte = (figD.Position(2)+figD.Position(4) + figD.Position(2))/2;
h.Label.String = '$F_5/5F_{\rm th}$';
h.Label.Interpreter = 'latex';
% set(ax,'YTick',0:0.2:0.8);
% set(ax,'XTick',[1,2,3],'XTickLabels',{'1','2','3'});
% set(ax,'YTick',0:0.2:0.8);
text(0.05,1.8,'(a)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);
set(ax,'FontWeight','bold');

%
ax = subplot(1,2,2);
blue = [0 82 177]/255;
red = [192 0 0]/255;
int = 1:0.01:(max(max(X))+0.01);
map = [1 1 1; blue];
map = interp1([0;1],map,linspace(0,1,length(int)));
map(1,:) = [];
contourf(gt/pi,gammat,X,int);
colormap(gca,map);
set(gca,'YTick',0:2);
set(gca,'XTick',0.1:0.4:0.9);
axis([0.1 0.9 0 2]);
xlabel('$g\tau_{\rm SA}/\pi$','interpreter','latex','FontSize',paper.labelSize );
ylabel('$\gamma\tau_{\rm SE}$','interpreter','latex','FontSize',paper.labelSize );
h=colorbar('eastoutside');
h.Ticks = [1 1.1];
h.Limits = [1 int(end)];
h.Label.String = '$F_5/5F_1$';
h.Label.Interpreter = 'latex';
set(ax,'FontSize',paper.fontSize);
set(ax,'FontWeight','bold');
text(0.15,1.8,'(b)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);

%
figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'TickLabelInterpreter','latex');
% figD.PaperPositionMode = 'auto';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
figD.Position(3) = paper.width*3;
% figD.Position(4) = paper.height;
% set(figD, 'color', 'none');
% %print('fig2.eps','-dpsc2');
print('figPartialSwap_Classical.pdf', '-dpdf', '-r600')
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);
