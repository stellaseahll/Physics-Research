clc; clear;clf;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure();

load('data1.mat');
loglog(nbar,Fg(:,end),'color',c(1,:),'linewidth',1)
hold on;
loglog(nbar,Fg(:,1)*np,'color',c(1,:),'linestyle',':','linewidth',1)
hold on;
loglog(nbar,Fe(:,end),'color',c(2,:),'linewidth',1)
hold on;
loglog(nbar,Fe(:,1)*np,'color',c(2,:),'linestyle',':','linewidth',1)

set(gca,'FontWeight','bold');
xlabel('\bf mean occ. num. $\bar{n} $','interpreter','latex','FontSize',8);
ylabel('\bf Fisher inf. $F$','interpreter','latex','FontSize',8);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1', '-dpdf', '-r600');
close(figD);

clc; clear;clf;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure();

load('data2.mat');
loglog(nbar,Fg(:,end),'color',c(1,:),'linewidth',1)
hold on;
loglog(nbar,Fg(:,1)*np,'color',c(1,:),'linestyle',':','linewidth',1)
hold on;
loglog(nbar,Fe(:,end),'color',c(2,:),'linewidth',1)
hold on;
loglog(nbar,Fe(:,1)*np,'color',c(2,:),'linestyle',':','linewidth',1)

set(gca,'FontWeight','bold');
xlabel('\bf mean occ. num. $\bar{n} $','interpreter','latex','FontSize',8);
ylabel('\bf Fisher inf. $F$','interpreter','latex','FontSize',8);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig2', '-dpdf', '-r600');
close(figD);

clc; clear;clf;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure();

load('data3.mat');
loglog(nbar,Fg(:,end),'color',c(1,:),'linewidth',1)
hold on;
loglog(nbar,Fg(:,1)*np,'color',c(1,:),'linestyle',':','linewidth',1)
hold on;
loglog(nbar,Fe(:,end),'color',c(2,:),'linewidth',1)
hold on;
loglog(nbar,Fe(:,1)*np,'color',c(2,:),'linestyle',':','linewidth',1)

set(gca,'FontWeight','bold');
xlabel('\bf mean occ. num. $\bar{n} $','interpreter','latex','FontSize',8);
ylabel('\bf Fisher inf. $F$','interpreter','latex','FontSize',8);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig3', '-dpdf', '-r600');
close(figD);
