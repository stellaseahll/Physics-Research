set(0,'defaulttextinterpreter','latex');
paper.width = 8;
paper.height = 2.5;
paper.fontSize = 12;
c1 = [0 82 177]/255;
c2 = [192 0 0]/255;
c3 = [84 130 53]/255;
c4 = [207 138 68]/255;
load('ohmic_xxx_partial_TC1_TH2_TW8.mat');
maxQC = max(max(QC));
xmax = 0.65;
ymax = 0.4;
interval =linspace(0,maxQC,11);
Ts = [TH TC TW];
nbar = 1./(exp(Es./Ts)-1);
keff = 1/6*kappa*sum(Es.*(2*nbar+1));
% Heating area
figD = figure();
[a b] = contour(g,keff,QC,[-0.1 0]);
x =a(1,2:105);
y =a(2,2:105);
hold on;
a(:,1) = [];
% a(:,end) = [];
heat = [x xmax xmax x(1); y y(end) 0  0];
red = c2;  
% Colormap
blue = c1;
darkblue = [48 86 102]/255;
k = [0;1];
map = [1 1 1; blue];
map = interp1(k,map,linspace(0,1,11));
map(1,:) = [];

colormap(map);
fill(heat(1,:),heat(2,:),red,'edgecolor','none')
% plot(a(1,:),a(2,:),'k','linewidth',2);
axis([0 xmax 0 ymax]);
contourf(g,keff,QC,interval,'linecolor','none');
set(gca,'XTick',0:0.2:0.6);
set(gca,'XTickLabel',{});
% xlabel('$g/ \Omega_{\mathrm {c}}$','interpreter','latex','FontSize',12);
ylabel('$\bar{\kappa}/\Omega_{\mathrm {c}}$','interpreter','latex','FontSize',12);
plot([0 0.397],[0 0.397],'color','k','linewidth',2);
kinterp = interp(keff,100);
for i = 1:length(g)
    QCtmp = interp(QC(:,i),100);
    [~, idx] = max(QCtmp);
    kmax(i) = kinterp(idx);
end
% [~ ,idx] = max(QCinterp);
plot(g(40:end-6),kmax(40:end-6),':','color','k','linewidth',2.2497);
%colorbar
h=colorbar('eastoutside');
% mitte = (figD.Position(2)+figD.Position(4) + figD.Position(2))/2;
% h.Label.String = '$\eta/\eta_{\mathrm{Otto}}$';
h.Label.Interpreter = 'latex';
set(get(h,'Label'),'string','$\mathcal{J}_{\mathrm {c}}^\infty/\hbar\Omega_{\mathrm {c}}^2$','interpreter','latex');
% h.Label.string = 0
h.Units = 'centimeters';
% set(h,'Ticks');
set(h,'FontWeight','bold');
set(h,'FontSize',12);
% axis([-3 1 -1 1]);

set(gca,'Units','centimeters');
set(gca,'FontWeight','bold');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);
print('QC', '-dpdf','-r600');
close(figD);