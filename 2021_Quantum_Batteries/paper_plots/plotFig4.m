%STEFAN NIMMRICHTER 29.4.2021
%Plot fig4 of paper

%load data
load('dataFig4.mat');
%[topt, cohEopt] = MEshorttime(0.5,200,100);

%%
figD = figure(); 
paper.width = 8;
paper.height = 3.5;
paper.fontSize = 7;
set(0,'defaulttextinterpreter','latex');
c = get(gca, 'ColorOrder');



fill([0 0 20 20 0],[0 1 1 0 0],0.9*[1 1 1],'LineStyle','none'); hold on;
thetaq = theta(1)*2;
% v1 = sin(thetaq).^2;
% v2 = sin(theta(1)).^2*(1-2*qopt(2));

y1 = linspace(0,numC*optCdE, length(cohE{1}(1,:)));
y2 = linspace(0,numC*optCdE, length(cohE{2}(1,:)));
y3 = topt*optCdE/optCAngle;

plot(topt,real(cohEopt(1,1:end-1))./y3,'linewidth',1.0,'color',[0 0 0]); hold on;
plot((1:numQ(2))*theta(2),cohE{2}(1,2:end)./y2(2:end),'Color',[132 47 64]/255,'MarkerfaceColor',[132 47 64]/255,'marker','^','linestyle','none','linewidth',1,'markersize',1.2); hold on;
plot((1:numQ(1))*theta(1),cohE{1}(1,2:end)./y1(2:end),'Color',[227 79 109]/255,'Markerfacecolor',[227 79 109]/255,'marker','s','linestyle','none','linewidth',1,'markersize',1.2); hold on;
plot((1:numQ(1))*thetaq,ones(length((1:numQ(1))*thetaq)),'color','k','marker','o','linestyle','none','linewidth',1,'markersize',1.2,'markerfacecolor',[0,0,0]);
hold off
set(gca,'FontSize',paper.fontSize);
% plot((0:numQ(2))*theta(2),v2*ones(length((0:numQ(2))*theta(2))),'color','k','marker','^','linestyle','none','linewidth',1.5,'markersize',2.5);
%set(gca, 'ColorOrder',c([3 1 2],:)); 
% plot((0:numQ(3))*theta(3),cohE{3}(1,:),'p','linewidth',1.5,'markersize',2); hold on;
% plot((0:numC)*optCAngle,(0:numC)*optCdE,'o-','linewidth',1.5,'markersize',2.5); hold on;
text(12,0.25,'classical processes','interpreter','latex','FontSize',paper.fontSize+2);
% plot((0:numQ(1))*theta(1),cohE{1}(1,:),'s-','linewidth',1.5,'markersize',2.5); hold on;
% plot((0:numQ(2))*theta(2),cohE{2}(1,:),'^-','linewidth',1.5,'markersize',2.5); hold on;
% plot(topt,real(cohEopt(1,1:end-1)),'k-','linewidth',1.5);hold on;
% plot((0:numQ(2))*theta(2),y2,'o');hold on;
% plot((0:numQ(1))*theta(1),y1);hold on;
% plot(topt,y3);

% for i = 1:3
%     plot((0:numQ(i))*theta(i),cohE{i}(1,:),'ro','linewidth',1.5,'markersize',2); hold on;
% end


axis([0 20 0 1.5]);
xlabel('charging time $k\theta$','interpreter','latex','FontSize',paper.fontSize+1);
%ylabel('power $\bar{n}_{c=1} (k)/\mathcal{P}_{c=0}^{\max} t$','interpreter','latex','FontSize',paper.fontSize);
ylabel('power $\mathcal{P}(k)/\mathcal{P}_{c=0}^{\max}$','interpreter','latex','FontSize',paper.fontSize+1);
set(gca,'Position',[0.10, 0.21, 0.88, 0.74]);
figD.Units = 'centimeters';
figD.PaperUnits = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
%set(figD.Children,'FontSize',paper.fontSize);
%figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig4', '-dpdf', '-r600','-loose');
%close(figD);
