paper.width = 5.5;
paper.height = 6.5;
paper.fontSize = 12;
set(0,'defaulttextinterpreter','latex')
load('fig1NEW.mat');
for i = 1:length(q)
    sEff2(i,2:tend+1) = sEff(i,2:tend+1)./(1:tend);
    tEff2(i,2:tend+1) = tEff(i,2:tend+1)./(1:tend);
    sEEff2(i,2:tend+1) = sEEff(i,2:tend+1)./(1:tend);
    tEEff2(i,2:tend+1) = tEEff(i,2:tend+1)./(1:tend);
end
% %% Contour plot
% figD = figure(); 
% contourf(1:tend,q,tEff2(:,2:end));
% set(gca,'FontWeight','bold');
% set(gca,'XTick',0:20:100);
% set(gca,'YTick',0:0.1:0.5);
% text(-22,0.48,'(a)','FontSize',paper.fontSize);
% h = colorbar;
% set(get(h,'label'),'string','efficiency $ \eta_{\rm erg}$','interpreter','latex');
% xlabel('$k$','interpreter','latex','FontSize',paper.fontSize);
% ylabel('$q$','interpreter','latex','FontSize',paper.fontSize);
% set(gca,'FontWeight','bold');
% % set(gca, 'XTickLabel', []);
% set(gca,'Units','centimeters');
% pos = get(gca,'Position');
% pos(3) = paper.width;
% pos(4) = paper.height;
% set(gca,'Position',pos);
% 
% %for saving
% set(figD.Children,'FontSize',paper.fontSize);
% print('fig1coh', '-dpdf', '-r600')
% close(figD);
% 
% %% Contour plot
% figD = figure(); 
% contourf(1:tend,q,sEff2(:,2:end));
% set(gca,'FontWeight','bold');
% set(gca,'XTick',0:20:100);
% set(gca,'YTick',0:0.1:0.5);
% text(-22,0.48,'(b)','FontSize',paper.fontSize);
% h = colorbar;
% set(get(h,'label'),'string','efficiency $ \eta_{\rm erg}$','interpreter','latex');
% xlabel('$k$','interpreter','latex','FontSize',paper.fontSize);
% ylabel('$q$','interpreter','latex','FontSize',paper.fontSize);
% set(gca,'FontWeight','bold');
% % set(gca, 'XTickLabel', []);
% set(gca,'Units','centimeters');
% pos = get(gca,'Position');
% pos(3) = paper.width;
% pos(4) = paper.height;
% set(gca,'Position',pos);
% 
% %for saving
% set(figD.Children,'FontSize',paper.fontSize);
% print('fig1incoh', '-dpdf', '-r600')
% close(figD);
load('acton10.mat', 'acton10');
acton10 = acton10(2:end,:);
acton10(end+1,:) = [1 1 1];
%% Contour plot
figD = figure(); 
sEff(end) = max(max(tEff));
contourf(0:tend,q,sEff);
colormap(gca,acton10);
% set(gca,'FontWeight','bold');
set(gca,'XTick',0:250:500);
set(gca,'YTick',[0:0.1:0.45 0.499]);
set(gca,'YTickLabel',0:0.1:0.5);
% text(600,-0.08,'(a)','FontSize',14);
text(-80,-0.14,'(a)','FontSize',paper.fontSize);
h = colorbar;
h.Ticks = 0:0.4:0.8;
h.TickLabels = 0:0.4:0.8;
h.Location = 'northoutside';
set(get(h,'label'),'string','$ \eta$','interpreter','latex');
xlabel('no. of qubits $k$','interpreter','latex','FontSize',paper.fontSize);
ylabel('population $q$','interpreter','latex','FontSize',paper.fontSize);
% set(gca,'FontWeight','bold');
% set(gca, 'XTickLabel', []);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1incoh', '-dpdf', '-r600')
close(figD);

figD = figure(); 
contourf(0:tend,q,tEff);
colormap(gca,acton10);
% set(gca,'FontWeight','bold');
set(gca,'XTick',0:250:500);
set(gca,'YTick',[0:0.1:0.45 0.499]);
set(gca,'YTickLabel',0:0.1:0.5);
% text(600,-0.08,'(a)','FontSize',14);
% text(-22,0.48,'(c)','FontSize',paper.fontSize);
text(-80,-0.14,'(b)','FontSize',paper.fontSize);
h = colorbar;
h.Ticks = 0:0.4:0.8;
h.TickLabels = 0:0.4:0.8;
h.Location = 'northoutside';
set(get(h,'label'),'string','$ \eta $','interpreter','latex');
xlabel('no. of qubits $k$','interpreter','latex','FontSize',paper.fontSize);
ylabel('population $q$','interpreter','latex','FontSize',paper.fontSize);
% set(gca,'FontWeight','bold');
% set(gca, 'XTickLabel', []);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1coh', '-dpdf', '-r600')
close(figD);

% paper.width = 10.5;
% paper.height = 6.5;
% paper.fontSize = 12;
%% Contour plot
figD = figure(); 
for i = 1:length(q)
[~,kActual(i)] = max(tEff(i,:));
end
theta = pi/2;
ptheta = sin(theta/2)^2;
c = sqrt(q.*(1-q));
N = 201;
kOpt = floor(N./(2*c*sqrt(ptheta*(1-ptheta))-ptheta*(2*q-1)));
% contourf(0:tend,q,tEff); hold on;
plot(kOpt,q,'k','linewidth',2);hold on;
plot(kActual,q,'k--','linewidth',2);
colormap(gca,acton10);
% set(gca,'FontWeight','bold');
axis([0 500 0 0.499]);
set(gca,'XTick',0:250:500);
set(gca,'YTick',[0:0.1:0.45 0.499]);
set(gca,'YTickLabel',0:0.1:0.5)
% text(600,-0.08,'(b)','FontSize',14);
h = colorbar;
h.Ticks = 0:0.4:0.8;
h.TickLabels = 0:0.4:0.8;
h.Location = 'northoutside';
set(get(h,'label'),'string',' $ \eta$','interpreter','latex');
xlabel('no. of qubits $k$','interpreter','latex','FontSize',paper.fontSize);
ylabel('population $q$','interpreter','latex','FontSize',paper.fontSize);
% set(gca,'FontWeight','bold');
% set(gca, 'XTickLabel', []);
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');

%for saving
set(figD.Children,'FontSize',paper.fontSize);
print('fig1curves', '-dpdf', '-r600')
close(figD);
