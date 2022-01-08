paper.width = 8;
paper.height = 9;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
% c = colormap(lines);
figD = figure(); 

load('1AncillaT2.mat');
p1th = nbar/(2*nbar+1);
FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1/(p1th-p1th^2);
partition = 0:0.2:2;
contourf(gammat,sin(gt),FFthRat,partition); hold on;
colormap(flipud(parula(length(partition)-1)));
axis([min(gammat) 2 min(sin(gt)) 1]);
set(gca,'XTick',[min(gammat) 1 2]);
set(gca,'XTickLabel',{0,1,2});
set(gca,'YTick',[min(sin(gt)) 0.5 1]);
set(gca,'YTickLabel',{0,0.5,1});
contour(gammat,sin(gt),FFthRat,0:1,'linewidth',1.5,'linestyle','--','color','k');
% shading interp

set(gca,'FontWeight','bold');
xlabel('\bf SE coupling $\kappa$ ','interpreter','latex','FontSize',8);
ylabel('\bf SA coupling $s$','interpreter','latex','FontSize',8);
text(1.8,0.1,'\bf (a)','interpreter','latex','FontSize',9)

h=colorbar('northoutside');
set(h, 'ylim', [0 4.035]);
set(h,'YTick',[0:0.5:2 4]);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1a.eps','-dpsc2');
close(figD);

paper.width = 8;
paper.height = 9;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
% c = colormap(lines);
% clf;
figD = figure(); 

load('1AncillaT20.mat');
p1th = nbar/(2*nbar+1);
FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1/(p1th-p1th^2);
FFthRat(end) = 4;
partition = 0:0.2:2;
contourf(gammat,sin(gt),FFthRat,partition); hold on;
colormap(flipud(parula(length(partition)-1)));
% contour(gammat,sin(gt),FFthRat,0:1,'linewidth',1.5,'linestyle','--','color','k');
axis([min(gammat) 2 min(sin(gt)) 1]);
set(gca,'XTick',[min(gammat) 1 2]);
set(gca,'XTickLabel',{0,1,2});
set(gca,'YTick',[min(sin(gt)) 0.5 1]);
set(gca,'YTickLabel',{0,0.5,1});
set(gca,'FontWeight','bold');
xlabel('\bf SE coupling $\kappa$ ','interpreter','latex','FontSize',8);
ylabel('\bf SA coupling $s$','interpreter','latex','FontSize',8);
text(1.8,0.1,'\bf (b)','interpreter','latex','FontSize',9)

h=colorbar('northoutside');
set(h, 'ylim', [0 4.035]);
set(h,'YTick',[0:0.5:2 4]);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1b.eps','-dpsc2');
close(figD);
