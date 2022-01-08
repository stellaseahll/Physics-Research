load('data5.mat');
str = 'tmp';
paper.width = 8;
paper.height = 5;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure(); 

plot(Dg,F(:,3),'.','color',c(3,:)); hold on;
plot(Dg,F(:,2),'.','color',c(2,:));
plot(Dg,F(:,1),'.','color',c(1,:));
set(gca,'FontWeight','bold');
axis tight
xlabel('\bf Geometric Discord $D_g$ ','interpreter','latex','FontSize',8);
ylabel('\bf $F_N$','interpreter','latex','FontSize',8);

figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print(strcat(str,'_Dg'), '-dpdf', '-r600')
close(figD);

figD = figure(); 
plot(I,F(:,3),'.','color',c(3,:)); hold on;
plot(I,F(:,2),'.','color',c(2,:));
plot(I,F(:,1),'.','color',c(1,:));
set(gca,'FontWeight','bold');
axis tight
xlabel('\bf Mutual Information $I$ ','interpreter','latex','FontSize',8);
ylabel('\bf $F_N$','interpreter','latex','FontSize',8);

figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print(strcat(str,'_MI'), '-dpdf', '-r600','-fillpage')
close(figD);
