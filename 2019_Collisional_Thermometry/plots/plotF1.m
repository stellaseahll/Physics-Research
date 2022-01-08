load('1AncillaT2_N500_gt_gammat.mat');
x1 = gammat;
y1 = gt;
%x1 = lamb;
%y1 = singt.^2;
z1 = FFthRat;
nbar1 = nbar;
%% 
load('1AncillaT2_N500_gt_gammat_plus.mat');
x2 = gammat;
y2 = gt;
%x2 = lamb;
%y2 = singt.^2;
z2 = FFthRat;
nbar2 = nbar;

%% plot

paper.width = 8;
paper.height = 5;
paper.fontSize = 8;
paper.labelSize = 10;
paper.xtickOff = 0.08; %manual offset of xTickLabels from axis
paper.xlabelOff = 0.15; %manual offset of xlabel from axis
set(0,'defaulttextinterpreter','latex');
% c = colormap(lines);
figD = figure(666); 
partition = 0:0.2:4;
colormap(flipud(hot(length(partition)-1)));

%ax1 = subplot(1,2,1,'align');
ax1 = subplot('position',[1/8 1/5 3.2/8 3.2/5]);
%fills w/o lines, then lines except at 1, then thick line at 1
[C1,H1] = contourf(x1,y1,z1,partition,'LineStyle','none'); hold on;
contour(x1,y1,z1,partition(partition~=1),'LineStyle','-','color','k');
set(ax1,'XTick',0:0.2:1);
%set(ax1,'XTickLabel',{''});
set(ax1,'YTick',0:0.2:1);
%set(ax1,'YTickLabel',{0,0.5,1});
contour(x1,y1,z1,[1 1],'linewidth',1.0,'linestyle',':','color','k');
% shading interp
hold off
caxis([partition(1) partition(end)]); %fix the color limits to partition limits!
axis([0 1 0 1]);

%manual xlabels to bring them closer to axis
% make a vector of vertical positions after the offset:
offset = repmat(ax1.YTick(1)-paper.xtickOff,1,numel(ax1.XTick));
% create new lables:
text(ax1.XTick,offset,ax1.XTickLabel,'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex');
ax1.XTickLabel=[]; %remove original ones

set(ax1,'FontSize',paper.fontSize);
%xl1 = xlabel('SE coupling $\lambda$','interpreter','latex','FontSize',paper.labelSize);
xl1 = xlabel('$\lambda$','interpreter','latex','FontSize',paper.labelSize);
xl1.Position(2) = -paper.xlabelOff;
%ylabel('SA coupling $s$','interpreter','latex','FontSize',paper.labelSize);
ylabel('$s$','interpreter','latex','FontSize',paper.labelSize);
text(0.75,0.15,'\bf (a)','interpreter','latex','FontSize',paper.labelSize)

%ax2 = subplot(1,2,2);
ax2 = subplot('position',[4.5/8 1/5 3.2/8 3.2/5]);
[C2,H2] = contourf(x2,y2,z2,partition); hold on;
contour(x2,y2,z2,partition(partition~=1),'LineStyle','-','color','k');
set(ax2,'XTick',0:0.2:1);
%set(ax2,'XTickLabel',{0,1,2});
set(ax2,'YTick',0:0.2:1);
set(ax2,'YTickLabel',[]);
contour(x2,y2,z2,[1 1],'linewidth',1.0,'linestyle',':','color','k');
% shading interp
hold off
caxis([partition(1) partition(end)]); %fix the color limits to partition limits!
axis([0 1 0 1]);

%manual xlabels to bring them closer to axis
% make a vector of vertical positions after the offset:
offset = repmat(ax2.YTick(1)-paper.xtickOff,1,numel(ax2.XTick));
% create new lables:
text(ax2.XTick,offset,ax2.XTickLabel,'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex');
ax2.XTickLabel=[]; %remove original ones

%set(ax2,'FontWeight','bold');
set(ax2,'FontSize',paper.fontSize);
%xl2 = xlabel('SE coupling $\lambda$','interpreter','latex','FontSize',paper.labelSize);
xl2 = xlabel('$\lambda$','interpreter','latex','FontSize',paper.labelSize);
xl2.Position(2) = -paper.xlabelOff;
ylabel([]);
text(0.75,0.15,'\bf (b)','interpreter','latex','FontSize',paper.labelSize)

h=colorbar('southoutside');
mitte = (ax2.Position(1)+ax2.Position(3) + ax1.Position(1))/2;
set(h,'Position',[mitte-2/8 4.4/5 4/8 0.2/5]);
set(h,'Ticks',0:0.2:4);
set(h,'TickLabels',{'0','','','','','1','','','','','2','','','','','3','','','','','4'});

%manual xlabels to bring them closer to axis
% make a vector of vertical positions after the offset:
offset2 = repmat(0.032,1,numel(h.XTick));
%FUCK: text uses units relative to axis object, cannot set to colorbar!
%--> convert normalized figure coordinates to normalized ax2 coord
Xfig2ax2=@(x)( (x-ax2.Position(1))/ax2.Position(3) );
Yfig2ax2=@(x)( (x-ax2.Position(2))/ax2.Position(4) );
% create new lables:
tt = text(Xfig2ax2(h.Position(1)+h.Position(3)*h.XTick/h.XTick(end)),Yfig2ax2(h.Position(2)+h.Position(4)+offset2),h.XTickLabel,...
    'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex','Units','normalized');
set(h,'TickLabels',[]); %remove original ones

%mark dotted line in colorbar at 1
L=annotation(figD,'line',(h.Position(1)+h.Position(3)/4)*[1 1],[h.Position(2) h.Position(2)+h.Position(4)],'LineStyle',':','LineWidth',1);

figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'TickLabelInterpreter','latex');
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
%print('fig1a.eps','-dpsc2');
%close(figD);
print('figF1', '-dpdf', '-r600')


