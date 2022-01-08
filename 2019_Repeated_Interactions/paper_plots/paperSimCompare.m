% Time evolutions for resonant exchange thermalization
% Paper plots
% Load paperSimCompare_getData to generate results

load('simResult2.mat')

%% Paper Plots a, b

paper.width = 8;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

c = colormap(lines);

figD = figure(); 
plot(gamma*t,real(obstS{1}(1,:)),'Color',c(1,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),real(obst{1}(1,1:30:end)),'Color',c(1,:),'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t,real(obstS{3}(1,:)),'Color',c(4,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),real(obst{3}(1,1:30:end)),'Color',c(4,:),'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t,real(obstS{2}(1,:)),'Color',c(2,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),real(obst{2}(1,1:30:end)),'Color',c(2,:),'Marker','^','MarkerSize',5,'LineStyle','none','LineWidth',1);

hold off
xlim([0,15]);
% ylim([-2,2]/omegaP);
ymin = -0.3;
ymax = 0;
% ylim([ymin ymax]);
set(gca,'FontWeight','bold')
% set(gca,'XTick',(0:5:20))
% set(gca,'XTickLabel',{1,10,100,1000})
% set(gca,'YTick',ymin:0.04:ymax)
xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% 
text(13,ymax-0.12*(ymax-ymin),'\bf (a)','interpreter','latex','FontSize',9)
% 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig4a', '-dpdf', '-r600')
% export_fig fig4a.pdf
close(figD)
% % 
% % 
paper.width = 8;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

c = colormap(lines);

figD = figure(); 

plot(gamma*t,abs(obstS{1}(3,:)),'Color',c(1,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),abs(obst{1}(3,1:30:end)),'Color',c(1,:),'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t,abs(obstS{2}(3,:)),'Color',c(2,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),abs(obst{2}(3,1:30:end)),'Color',c(2,:),'Marker','^','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t,abs(obstS{3}(3,:)),'Color',c(4,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),abs(obst{3}(3,1:30:end)),'Color',c(4,:),'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);


hold off
% xlim(gamma*[20,4e5]);
% ylim([-2,2]/omegaP);
ymin = 0;
ymax = 0.5;
xlim([0,15]);
ylim([ymin ymax]);
set(gca,'FontWeight','bold')
% set(gca,'XTick',(0:2:10))
% set(gca,'XTickLabel',{1,10,100,1000})
% set(gca,'YTick',ymin:0.2:ymax)
xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
ylabel('\bf coherence $|\langle \hat{J}_+\rangle|$','interpreter','latex','FontSize',8)
% % 
text(13,ymax-0.12*(ymax-ymin),'\bf (b)','interpreter','latex','FontSize',9)
% % 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig4b', '-dpdf', '-r600')
% export_fig fig4a.pdf
close(figD)
% 
paper.width = 8;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

c = colormap(lines);

figD = figure(); 
plot(gamma*t,real(obstS{1}(1,:)),'Color',c(1,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),real(obstEik{1}(1,1:30:end)),'Color',c(2,:),'Marker','^','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t(1:30:end),real(obstUint{1}(1,1:30:end)),'Color',c(4,:),'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);

hold off
xlim([0,15]);
% ylim([-2,2]/omegaP);
ymin = -0.3;
ymax = 0;
% ylim([ymin ymax]);
set(gca,'FontWeight','bold')
% set(gca,'XTick',(0:5:20))
% set(gca,'XTickLabel',{1,10,100,1000})
% set(gca,'YTick',ymin:0.04:ymax)
xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% 
text(13,ymax-0.12*(ymax-ymin),'\bf (a)','interpreter','latex','FontSize',9)
% 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig5a', '-dpdf', '-r600')
% export_fig fig4a.pdf
close(figD)
% % 
paper.width = 8;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

c = colormap(lines);

figD = figure(); 
plot(gamma*t,real(obstS{2}(1,:)),'Color',c(1,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),real(obstEik{2}(1,1:30:end)),'Color',c(2,:),'Marker','^','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t(1:30:end),real(obstUint{2}(1,1:30:end)),'Color',c(4,:),'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);

hold off
xlim([0,15]);
% ylim([-2,2]/omegaP);
ymin = -0.3;
ymax = 0;
% ylim([ymin ymax]);
set(gca,'FontWeight','bold')
% set(gca,'XTick',(0:5:20))
% set(gca,'XTickLabel',{1,10,100,1000})
% set(gca,'YTick',ymin:0.04:ymax)
xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% 
text(13,ymax-0.12*(ymax-ymin),'\bf (b)','interpreter','latex','FontSize',9)
% 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig5b', '-dpdf', '-r600')
% export_fig fig4a.pdf
close(figD)
% % 

% paper.width = 8;
% paper.height = 4;
% paper.fontSize = 8;
% set(0,'defaulttextinterpreter','latex')
% 
% c = colormap(lines);
% 
% figD = figure(); 
% plot(gamma*t,real(obstS{3}(1,:)),'Color',c(1,:),'LineWidth',1); hold on;
% plot(gamma*t(1:30:end),real(obstEik{3}(1,1:30:end)),'Color',c(2,:),'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1);
% plot(gamma*t(1:30:end),real(obstUint{3}(1,1:30:end)),'Color',c(4,:),'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1);
% 
% hold off
% xlim([0,15]);
% % ylim([-2,2]/omegaP);
% ymin = -0.3;
% ymax = 0;
% % ylim([ymin ymax]);
% set(gca,'FontWeight','bold')
% % set(gca,'XTick',(0:5:20))
% % set(gca,'XTickLabel',{1,10,100,1000})
% % set(gca,'YTick',ymin:0.04:ymax)
% xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
% ylabel('\bf mean spin $\langle\hat{J}_z\rangle$','interpreter','latex','FontSize',8)
% % 
% text(13,ymax-0.12*(ymax-ymin),'\bf (g)','interpreter','latex','FontSize',9)
% % 
% %for saving
% figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'FontSize',paper.fontSize)
% figD.PaperPositionMode = 'auto';
% posD=get(figD,'Position');
% set(figD,'PaperSize',[posD(3), posD(4)])
% set(figD, 'color', 'none');
% print('fig4g', '-dpdf', '-r600')
% % export_fig fig4a.pdf
% close(figD)

%%
paper.width = 8;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

c = colormap(lines);

figD = figure(); 
plot(gamma*t,abs(obstS{1}(3,:)),'Color',c(1,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),abs(obstEik{1}(3,1:30:end)),'Color',c(2,:),'Marker','^','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t(1:30:end),abs(obstUint{1}(3,1:30:end)),'Color',c(4,:),'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);

hold off
ymin = 0;
ymax = 0.5;
xlim([0,15])
ylim([ymin ymax]);
set(gca,'FontWeight','bold')
% set(gca,'XTick',(0:2:10))
% set(gca,'XTickLabel',{1,10,100,1000})
% set(gca,'YTick',ymin:0.2:ymax)
xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
ylabel('\bf coherence $|\langle \hat{J}_+\rangle|$','interpreter','latex','FontSize',8)
% % 
text(13,ymax-0.12*(ymax-ymin),'\bf (c)','interpreter','latex','FontSize',9)
% 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig5c', '-dpdf', '-r600')



%%

paper.width = 8;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex')

c = colormap(lines);

figD = figure(); 
plot(gamma*t,abs(obstS{2}(3,:)),'Color',c(1,:),'LineWidth',1); hold on;
plot(gamma*t(1:30:end),abs(obstEik{2}(3,1:30:end)),'Color',c(2,:),'Marker','^','MarkerSize',5,'LineStyle','none','LineWidth',1);
plot(gamma*t(1:30:end),abs(obstUint{2}(3,1:30:end)),'Color',c(4,:),'Marker','*','MarkerSize',5,'LineStyle','none','LineWidth',1);

hold off
ymin = 0;
ymax = 0.5;
xlim([0,15])
ylim([ymin ymax]);
set(gca,'FontWeight','bold')
% set(gca,'XTick',(0:2:10))
% set(gca,'XTickLabel',{1,10,100,1000})
% set(gca,'YTick',ymin:0.2:ymax)
xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
ylabel('\bf coherence $|\langle \hat{J}_+\rangle|$','interpreter','latex','FontSize',8)
% % 
text(13,ymax-0.12*(ymax-ymin),'\bf (d)','interpreter','latex','FontSize',9)
% 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize)
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print('fig5d', '-dpdf', '-r600')

%%
% 
% paper.width = 8;
% paper.height = 4;
% paper.fontSize = 8;
% set(0,'defaulttextinterpreter','latex')
% 
% c = colormap(lines);
% 
% figD = figure(); 
% plot(gamma*t,abs(obstS{3}(3,:)),'Color',c(1,:),'LineWidth',1); hold on;
% plot(gamma*t(1:30:end),abs(obstEik{3}(3,1:30:end)),'Color',c(2,:),'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1);
% plot(gamma*t(1:30:end),abs(obstUint{3}(3,1:30:end)),'Color',c(4,:),'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1);
% 
% hold off
% ymin = 0;
% ymax = 0.5;
% xlim([0,15])
% ylim([ymin ymax]);
% set(gca,'FontWeight','bold')
% % set(gca,'XTick',(0:2:10))
% % set(gca,'XTickLabel',{1,10,100,1000})
% % set(gca,'YTick',ymin:0.2:ymax)
% xlabel('\bf time $\gamma t$','interpreter','latex','FontSize',8)
% ylabel('\bf coherence $|\langle \hat{J}_+\rangle|$','interpreter','latex','FontSize',8)
% % % 
% text(13,ymax-0.12*(ymax-ymin),'\bf (d)','interpreter','latex','FontSize',9)
% % 
% %for saving
% figD.Units = 'centimeters';
% figD.Position(3) = paper.width;
% figD.Position(4) = paper.height;
% set(figD.Children,'FontSize',paper.fontSize)
% figD.PaperPositionMode = 'auto';
% posD=get(figD,'Position');
% set(figD,'PaperSize',[posD(3), posD(4)])
% set(figD, 'color', 'none');
% print('fig5d', '-dpdf', '-r600')
