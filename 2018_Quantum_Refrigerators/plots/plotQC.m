clf;clear;
load('ohmic_xxx_partial_TC1_TH2_TW8.mat');
QC = QC{1};
QC(1) = 0.5;
xmax = 0.65;
ymax = 0.4;
interval = 0:0.00025:0.00325;
Ts = [TH TC TW];
nbar = 1./(exp(Es./Ts)-1);
keff = 1/6*kappa*sum(Es.*(2*nbar+1));
% Heating area
[a b] = contour(g,keff,QC,[-0.1 0]);
% hold on;
a(:,1) = [];
heat = [a(1,end:-1:1) xmax xmax a(1,1); a(2,end:-1:1) 0 ymax ymax];
red = [202 76 109 ]/255;
% Colormap
blue = [168 230 248]/255;
darkblue = [4 65 81]/255;
k = [0; 0.4;1];
map = [1 1 1; blue; darkblue];
map = interp1(k,map,linspace(0,1,length(interval)-1));

% color = zeros(floor(length(interval)/2)+1,4);
% color(:,1) = 0.53;
% color(:,2) = 0.16;
% color(:,3) = 0;
% color(:,4) = 5.5/100*[0:floor(length(interval)/2)]';
% map = cmyk2rgb(color);
% map2 = [1 1 1; map(1,:)];
% k = [0;1];
% map2 = interp1(k,map2,linspace(0,1,2+floor(length(interval)/2)));
% map = [map2(3:end-1,:) ;map];
% k = interval-min(interval);
% k = interval/max(interval);
% map = interp1(k,map);
contourf(g,keff,QC,interval,'linewidth',2);
hold on;
colormap(map);
fill(heat(1,:),heat(2,:),red,'edgecolor','none')
plot(a(1,:),a(2,:),'k','linewidth',2);
blue = [46 0 202]/255;
red = [227 0 115]/255;
plot([0 ymax],[0 ymax],'-.','color',blue,'linewidth',3);
kinterp = interp(keff,100);
for i = 1:length(g)
    QCtmp = interp(QC(:,i),100);
    [~, idx] = max(QCtmp);
    kmax(i) = kinterp(idx);
end
% [~ ,idx] = max(QCinterp);
plot(g(40:end-6),kmax(40:end-6),':','color',red,'linewidth',3);
axis([0 xmax 0 ymax]);
set(gca,'FontSize', 18);
set(gca,'XTick',0:0.2:0.6);
set(gca,'YTick',0:0.1:0.4);
xlabel('T_W-T_C');
ylabel('T_H-T_C');
h = colorbar;
set(h, 'ylim', [0 0.00325])
set(h,'YTick',0:0.001:0.003);
print('test.eps','-dpsc2');