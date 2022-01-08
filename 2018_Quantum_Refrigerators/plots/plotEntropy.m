clf;
load('ohmic_xxx_partial_TC1_TH2_TW8.mat');
% QC = QC{1};
% QC(1) = 0.5;
Smax = max(max(S));
S(end,end) = 0.12;
Smin = min(min(S));
interval = 0:0.02:Smax+0.01;
% Ts = [TH TC TW];
Es = [5 1 4];
Ts = [2 1 8];
nbar = 1./(exp(Es./Ts)-1);
chi = logspace(-3,-1,100);
keff = 1/6*chi*sum(Es.*(2*nbar+1));
% Heating area
[a b] = contour(g,keff,S,interval);
% hold on;
% a(:,1) = [];
% red = [251,79,79]/255;
% % Colormap
% Colormap
blue = [108,192,229]/255;
darkblue = [48 86 102]/255;
k = [0;1];
map = [1 1 1; darkblue];
map = interp1(k,map,linspace(0,1,length(interval)-1));
% map = [map2(3:end-1,:) ;map];
% k = interval-min(interval);
% k = interval/max(interval);
contourf(g,keff,S,interval,'linewidth',2);
hold on;
colormap(map);
% fill(heat(1,:),heat(2,:),red,'edgecolor','none')
% plot(a(1,:),a(2,:),'k','linewidth',2);
% blue = [46 0 202]/255;
% red = [227 0 115]/255;
% plot([0 ymax],[0 ymax],'-.','color',blue,'linewidth',3);
% kinterp = interp(keff,100);
% for i = 1:length(g)
%     QCtmp = interp(QC(:,i),100);
%     [~, idx] = max(QCtmp);
%     kmax(i) = kinterp(idx);
% end
% [~ ,idx] = max(QCinterp);
% plot(g(40:end-6),kmax(40:end-6),':','color',red,'linewidth',3);
% axis([0 xmax 0 ymax]);
set(gca,'FontSize', 18);
set(gca,'XTick',0:0.2:1);
set(gca,'YTick',0:0.1:0.4);
% xlabel('T_W-T_C');
% ylabel('T_H-T_C');
h = colorbar;
set(h, 'ylim', [0 Smax])
set(h,'YTick',0:0.02:0.12);
print('test.eps','-dpsc2');