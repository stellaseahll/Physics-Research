paper.width = 8;
paper.height = 2.5;
paper.fontSize = 12;
set(0,'defaulttextinterpreter','latex');
load('thesisIncoherentData_wp1_ws100_g2.50_kh0.0010_nh1.0_nc0.1.mat');
work = -real(W)/ws/kh;
eff = real(-W)./real(Qh);

c1 = [0 82 177]/255;
c2 = [192 0 0]/255;
c3 = [84 130 53]/255;
c4 = [207 138 68]/255;
red = c2;
blue = c1;
map = [1 1 1; red];
k = [0; 1];
map = interp1(k,map,linspace(0,1,13));
% map(1,:) = [];

figD = figure;
[a b] = contour(log10(km),log10(kc),work,[-0.1 0]);
x =a(1,2:end);
y =a(2,2:end);
negwork = [x -1 -1 -3;y y(end) -3 -3];
red = c2;  
fill(negwork(1,:),negwork(2,:),blue,'edgecolor','none');
hold on;
contourf(log10(km),log10(kc),work,0:0.05:0.7,'linecolor','none');
colormap(map);
ylabel('$\kappa_{\rm c}/\omega$','interpreter','latex','FontSize',12);
% xlabel('$\gamma/\omega$','interpreter','latex','FontSize',12);
set(gca,'XTick',-3:-1,'XTickLabels',{});
set(gca,'YTick',-3:-1,'YTickLabels',{'0.001','0.01','0.1'});
set(gca,'FontWeight','bold');
% text(0.1,0.92,'(a)','interpreter','latex','FontSize',12);
% set(ax,'FontSize',paper.fontSize);


h=colorbar('eastoutside');
% axis([-3 1 -1 1]);
% mitte = (figD.Position(2)+figD.Position(4) + figD.Position(2))/2;
% h.Label.String = '$\eta/\eta_{\mathrm{Otto}}$';
h.Label.Interpreter = 'latex';
set(get(h,'Label'),'string','$\mathcal{W}/\hbar\Omega\kappa_{\mathrm{h}}$','interpreter','latex');
% h.Label.string = 0
h.Units = 'centimeters';
set(h,'Ticks',0:0.2:0.8);
set(h,'FontWeight','bold');
set(h,'FontSize',12);


set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);

print('incoherentWorkContour', '-dpdf', '-r600')
close(figD);
figD = figure();
darkred = [117 24 16]/255;
map = [1 1 1; red];
k = [0; 1];
map = interp1(k,map,linspace(0,1,9));
fill(negwork(1,:),negwork(2,:),blue,'edgecolor','none');
hold on;
contourf(log10(km),log10(kc),eff,0:0.1:0.8,'linecolor','none');
colormap(map);
ylabel('$\kappa_{\rm c}/\omega$','interpreter','latex','FontSize',12);
xlabel('$\gamma/\omega$','interpreter','latex','FontSize',12);
set(gca,'XTick',-3:-1,'XTickLabels',{'0.001','0.01','0.1'});
set(gca,'YTick',-3:-1,'YTickLabels',{'0.001','0.01','0.1'});
set(gca,'FontWeight','bold');
% text(0.1,0.92,'(a)','interpreter','latex','FontSize',12);
% set(ax,'FontSize',paper.fontSize);


h=colorbar('eastoutside');
caxis([0 0.9]);
% axis([-3 1 -1 1]);
% mitte = (figD.Position(2)+figD.Position(4) + figD.Position(2))/2;
% h.Label.String = '$\eta/\eta_{\mathrm{Otto}}$';
h.Label.Interpreter = 'latex';
set(get(h,'Label'),'string','$\eta$','interpreter','latex');
% h.Label.string = 0
h.Units = 'centimeters';
set(h,'Ticks',[0.00001 0.4  0.8],'TickLabels',{0  0.4  0.8});
set(h,'FontWeight','bold');
set(h,'FontSize',12);


set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);

print('incoherentEffContour', '-dpdf', '-r600')
close(figD);