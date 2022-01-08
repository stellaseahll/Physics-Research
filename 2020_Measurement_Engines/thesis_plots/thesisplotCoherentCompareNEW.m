paper.width = 8;
paper.height = 2.5;
paper.fontSize = 12;
set(0,'defaulttextinterpreter','latex')

%here saved as arrays for a fixed gt=pi/100
% load('CoherentData_Global_Step_wp1_ws100_g2.50_G0.1000_kh0.0010_kc0.1000_nc0.10.mat'); d1 = delta/(g*g/wp)*2; w1 = real(-W)/kh/ws; n1 = real(-W)./real(Qh); 
% load('CoherentData_Global_wp1_ws100_g2.50_G0.1000_kh0.0010_kc0.1000_nc0.10.mat'); d2 = delta/(g*g/wp)*2; w2 = real(-W)/kh/ws; n2 = real(-W)./real(Qh); 
load('CoherentDataNEW3_Global_Step_wp1_ws100_g3.54_G0.1000_kh0.0010_kc0.1000_nc1.00.mat'); d1 = delta/(g*g/wp)*2; w1 = real(-W)/kh/ws; n1 = real(-W)./real(Qh); 
load('CoherentDataNEW3_Global_wp1_ws100_g3.54_G0.1000_kh0.0010_kc0.1000_nc1.00.mat'); d2 = delta/(g*g/wp)*2; w2 = real(-W)/kh/ws; n2 = real(-W)./real(Qh); 
% load('Projector_nc1.000.mat'); d3 = delta/(g*g/wp); w3 = real(-W)/kh/ws; n3 = real(-W)./real(Qh); 
% x0 = sqrt(2)*g;
% Wmax = G*(ws-2*wp*x0^2)/kh/ws*nh/(2*nh+1);
n1(w1<0) = NaN;
n2(w2<-0.1) = NaN;
% n3(w3<0) = 0;
%% plot wide fig (small g) 

figD = figure();
c1 = [52 79 168]/255;
c2 = [193 24 0]/255;
c3 = [84 130 53]/255;
c4 = [207 138 68]/255;
c = [c1;c2;c3;c4];


% a = 0.9;
% fill([0 0 XEnd XEnd 0], [-0.1 0.8 0.8 -0.1 -0.1],[a a a],'linestyle','none');
plot(d1,w1,'Color',c1,'LineWidth',2);hold on;
% plot(d2,w2,'Color',c(2,:),'LineWidth',1);
plot(d2,w2,'Color',c2,'LineWidth',2);
% plot([-8 8],[Wmax Wmax],'Color',c(1,:),'LineWidth',1,'LineStyle',':');
% plot(X,w3,'Color',c(4,:),'LineWidth',1);
% plot([0 X(end)],[wth1 wth1],':','Color',c(1,:),'LineWidth',1);
% plot([0 X(end)],[wth2 wth2],':','Color',c(2,:),'LineWidth',1);
% plot([0 X(end)],[wth3 wth3],':','Color',c(4,:),'LineWidth',1);
xlabel('\bf Detuning $\Delta/\omega x_0^2$','interpreter','latex','FontSize',12);
ylabel('\bf Power $\mathcal{W}/\hbar\Omega\kappa_{\rm h}$','interpreter','latex','FontSize',12);
xlim([-8 8]);
ylim([0 0.4]);
set(gca,'XTick',-8:4:8);
set(gca,'YTick',0:0.2:0.4);
set(gca,'FontWeight','bold');
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);

print('coherentWork', '-dpdf', '-r600')
close(figD);
% set(ax,'XTick',[1,2,3],'XTickLabels',{'1','2','3'});
% set(ax,'YTick',0:0.2:0.8);
% text(-7.0,0.34,'(a)','interpreter','latex','FontSize',12);
% set(ax,'FontSize',paper.fontSize);
figD = figure;
plot(d1,n1,'Color',c1,'LineWidth',2);hold on;
% plot(d2,w2,'Color',c(2,:),'LineWidth',1);
plot(d2,n2,'Color',c2,'LineWidth',2);
xlabel('\bf Detuning $\Delta/\omega x_0^2$','interpreter','latex','FontSize',12);
ylabel('\bf Efficiency $\eta$','interpreter','latex','FontSize',12);

xlim([-8 8]);
ylim([0 1]);
set(gca,'XTick',-8:4:8);
set(gca,'YTick',0:0.5:1);
set(gca,'FontWeight','bold');
% 
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);

print('coherentEff', '-dpdf', '-r600')

% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);