clear;
gt = pi/2;
gammat = logspace(-2,1,100);
nbar = logspace(-2,2,100);
dnbar = 1e-5*nbar;
Fth = 1./(1+2*nbar).^2./(1+nbar)./nbar;
Fth = Fth'*ones(1,100);
rho = [0 0; 0 1];
for i = 1:length(nbar)
    fprintf('(%d)\n',i);
    for j = 1:length(gammat)
        s = spSinglePassXYNEW(gt,gammat(j),nbar(i),dnbar(i),rho,2,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
    end
end
paper.width = 8;
paper.height = 9;
paper.fontSize = 14;
set(0,'defaulttextinterpreter','latex');
% c = colormap(lines);
figD = figure(); 


partition =-2:0.5:3.5;
c0 = [255,255,255; 255,237,140; 254,178,76; 240,59,32; 0,0,0]/255;
cm = interp1([0;1/5;1/1.8;3/4;1.2],c0,linspace(0,1,length(partition)-1));
colormap(cm);
FFth = log10(F1./Fth);
FFth(FFth==min(min(FFth))) = -2;
FFth(FFth==max(max(FFth))) = 3.5;
contourf(log10(gammat),log10(nbar),FFth,partition); hold on;
% colormap(flipud(parula(length(partition)-1)));
axis([min(log10(gammat)) max(log10(gammat)) min(log10(nbar)) max(log10(nbar))]);
set(gca,'XTick',-2:1);
set(gca,'XTickLabel',{-2,-1,0,1});
set(gca,'YTick',-2:2);
% set(gca,'YTickLabel',{0,0.5,1});
contour(log10(gammat),log10(nbar),FFth,[0 10],'linewidth',1.2,'linestyle','--','color','k');
% shading interp

set(gca,'FontWeight','bold');
xlabel('$\log(\gamma\tau_{\mathrm{SE}})$ ','interpreter','latex','FontSize',14);
ylabel('$\log\bar{n}$','interpreter','latex','FontSize',14);

h=colorbar('northoutside');
set(h, 'ylim', [-2 3.5]);
set(h,'YTick',[-2:2:3.5]);
h.Label.String = '$\log (F_1/F_{\mathrm{th}})$';
h.Label.FontSize = 14;
h.Label.Interpreter = 'latex';
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1aex.eps','-dpsc2');
close(figD);