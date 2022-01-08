%% load data
panels = {'(a)','(b)'};

% load('8AncillaT2_g_e_p_tol6.mat');
% Fg = Fs{1}(1,:);
% Fe = Fs{2}(1,:);
% Fp = Fs{3}(1,:);

%here saved as arrays for a fixed gt=pi/100
load('CoherentData_Global_Step_wp1_ws100_g2.50_G0.1000_kh0.0010_kc0.1000_nc0.10.mat'); d1 = delta/(g*g/wp)*2; w1 = real(-W)/kh/ws; n1 = real(-W)./real(Qh); 
load('CoherentData_Global_Step_wp1_ws100_g2.50_G0.1000_kh0.0010_kc0.1000_nc1.mat'); d2 = delta/(g*g/wp)*2; w2 = real(-W)/kh/ws; n2 = real(-W)./real(Qh); 
% load('Projector_nc1.000.mat'); d3 = delta/(g*g/wp); w3 = real(-W)/kh/ws; n3 = real(-W)./real(Qh); 
n1(w1<0) = 0;
n2(w2<0) = 0;
% n3(w3<0) = 0;
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
plot(d1,w1,'Color',c(1,:),'LineWidth',1);hold on;
% plot(d2,w2,'Color',c(2,:),'LineWidth',1);
plot(d2,w2,'Color',c(2,:),'LineWidth',1);
% plot(X,w3,'Color',c(4,:),'LineWidth',1);
% plot([0 X(end)],[wth1 wth1],':','Color',c(1,:),'LineWidth',1);
% plot([0 X(end)],[wth2 wth2],':','Color',c(2,:),'LineWidth',1);
% plot([0 X(end)],[wth3 wth3],':','Color',c(4,:),'LineWidth',1);
xlabel('Detuning $\Delta/x_0^2\omega$','interpreter','latex','FontSize',paper.labelSize);
ylabel('Work power $\dot{W}/\hbar\Omega\kappa_h$','interpreter','latex','FontSize',paper.labelSize);
xlim([-8 8]);
ylim([-0.05 0.8]);
set(ax,'XTick',-8:2:8);
set(ax,'YTick',0:0.2:0.8);
% set(ax,'XTick',[1,2,3],'XTickLabels',{'1','2','3'});
% set(ax,'YTick',0:0.2:0.8);
text(-3.6,0.65,'(a)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);

ax = subplot(1,2,2);
plot(d1,n1,'Color',c(1,:),'LineWidth',1);hold on;
% plot(d2,n2,'Color',c(2,:),'LineWidth',1);
plot(d2,n2,'Color',c(2,:),'LineWidth',1);
xlabel('Detuning $\Delta/x_0^2\omega$','interpreter','latex','FontSize',paper.labelSize);
ylabel('Efficiency $\eta$','interpreter','latex','FontSize',paper.labelSize);

xlim([-8 8]);
ylim([0 1]);
set(ax,'XTick',-8:2:8);
set(ax,'YTick',0:0.2:1);
text(-3.6,0.82,'(b)','interpreter','latex','FontSize',paper.labelSize);
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
print('fig2.pdf', '-dpdf', '-r600')
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);

