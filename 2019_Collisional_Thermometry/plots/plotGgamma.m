clc; clear;
paper.width = 7;
paper.height = 4;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure(); 

% close(figure);

%% Plot contour for max Fisher and alpha
theta = 0.5;
nbar = 1;
np = 6;
isSteady = 1;
filename = sprintf('theta%.1f,nbar%.3f,np%d,isSteady%d.mat',theta,nbar,np,isSteady);
load(filename);
for i = 1:length(G)
    for j = 1:length(gammat)
        Fx(i,j) = F(np,i,j);
    end
end

contourf(1-exp(-gammat*(2*nbar+1)),sin(G(21:40)*pi).^2,Fx(21:40,:));


d= colorbar;
set(gca,'FontWeight','bold')
set(gca,'XTick',[min(1-exp(-gammat*(2*nbar+1))),0.5,max(1-exp(-gammat*(2*nbar+1)))]);
set(gca,'XTickLabel',{0,0.5,1});
set(gca,'YTick',[min(sin(G(21:40)*pi).^2),0.5,max(sin(G(21:40)*pi).^2)]);
set(gca,'YTickLabel',{0,0.5,1});
xlabel('\bf bath weight $\Lambda$ ','interpreter','latex','FontSize',8);
ylabel('\bf interaction weight $K$','interpreter','latex','FontSize',8);
d.Label.Interpreter = 'latex';
str = sprintf('$F_%d$',np);
d.Label.String = str;
d.Label.FontSize = 8;
d.Ticks = 0:0.04:1; 
%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig3a', '-dpdf', '-r600')
close(figD);

%%

figD = figure(); 
contourf(1-exp(-gammat*(2*nbar+1)),sin(G(21:40)*pi).^2,real(alpha(21:40,:)));
d= colorbar;
set(gca,'FontWeight','bold')
set(gca,'XTick',[min(1-exp(-gammat*(2*nbar+1))),0.5,max(1-exp(-gammat*(2*nbar+1)))]);
set(gca,'XTickLabel',{0,0.5,1});
set(gca,'YTick',[min(sin(G(21:40)*pi).^2),0.5,max(sin(G(21:40)*pi).^2)]);
set(gca,'YTickLabel',{0,0.5,1});
xlabel('\bf bath weight $\Lambda$ ','interpreter','latex','FontSize',8);
ylabel('\bf interaction weight $K$','interpreter','latex','FontSize',8);
d.Label.Interpreter = 'latex';
d.Label.String = '$\alpha$';
d.Label.FontSize = 8;
d.Ticks = 0:2:10; 

%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
% print('fig3b', '-dpdf', '-r600')
close(figD);