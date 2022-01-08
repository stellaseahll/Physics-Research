%STEFAN NIMMRICHTER 29.4.2021
%Plot fig5 of paper

load('fig1NEW.mat');
for i = 1:length(q)
    sEff2(i,2:tend+1) = sEff(i,2:tend+1)./(1:tend);
    tEff2(i,2:tend+1) = tEff(i,2:tend+1)./(1:tend);
    sEEff2(i,2:tend+1) = sEEff(i,2:tend+1)./(1:tend);
    tEEff2(i,2:tend+1) = tEEff(i,2:tend+1)./(1:tend);
end

load('acton10.mat', 'acton10');
acton10 = acton10(2:end,:);
acton10(end+1,:) = [1 1 1];

for i = 1:length(q)
[~,kActual(i)] = max(tEff(i,:));
end
theta = pi/2;
ptheta = sin(theta/2)^2;
c = sqrt(q.*(1-q));
N = 201;
kOpt = floor(N./(2*c*sqrt(ptheta*(1-ptheta))-ptheta*(2*q-1)));

%% plot

paper.width = 8;
paper.height = 3.7;
paper.bottomgap = 0.7;
paper.fontSize = 7;
paper.labelSize = 8;
paper.xtickOff = 0.065; %manual offset of xTickLabels from axis
paper.xlabelOff = 0.09; %manual offset of xlabel from axis
set(0,'defaulttextinterpreter','latex');

figD = figure(666); 

ax1pos = [0.98/paper.width paper.bottomgap/paper.height 2.9/paper.width 2.9/paper.height];
ax2pos = ax1pos; ax2pos(1) = 4.14/paper.width;
hpos = [7.3/paper.width paper.bottomgap/paper.height 0.2/paper.width 2.1/paper.height];

partition = 0:0.1:0.8;
colormap(acton10(2:end,:));

sEff2 = sEff;
sEff2(end) = max(max(tEff));
ax1 = subplot('position',ax1pos);
[C1,H1] = contourf(0:tend,q,sEff2,partition,'ShowText','off'); hold on;
%caxis([0 0.4]); %fix the color limits to partition limits!
set(ax1,'XTick',0:100:500);
%set(ax1,'XTickLabel',{'0','','200','','400',''});
set(ax1,'XTickLabel',{'','100','','300','','500'});
set(ax1,'YTick',0:0.1:0.5);
%set(ax1,'YTickLabel',{'0','0.1','0.2','0.3','0.4',''});
hold off
axis([0 500 0 0.5]);

%manual xlabels to bring them closer to axis
% make a vector of vertical positions after the offset:
%offset = repmat(ax1.YTick(1)-paper.xtickOff,1,numel(ax1.XTick));
% create new lables:
%text(ax1.XTick,offset,ax1.XTickLabel,'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex');
%ax1.XTickLabel=[]; %remove original ones

set(ax1,'FontSize',paper.fontSize);
%xl1 = xlabel('SE coupling $\lambda$','interpreter','latex','FontSize',paper.labelSize);
xl1 = xlabel('no.~of qubits $k$','interpreter','latex','FontSize',paper.labelSize);
%xl1.Position(2) = -paper.xlabelOff;
%ylabel('SA coupling $s$','interpreter','latex','FontSize',paper.labelSize);
ylabel('population $q$','interpreter','latex','FontSize',paper.labelSize);
ta = text(-168,0.48,'(a)','interpreter','latex','FontSize',paper.labelSize+1);

%ax2 = subplot(1,2,2);
ax2 = subplot('position',ax2pos);
[C2,H2] = contourf(0:tend,q,tEff,partition(2:end),'ShowText','off'); hold on;
%contour(x2,y2,z2,partition(partition~=1),'LineStyle','-','color','k','linewidth',0.25);
set(ax2,'XTick',0:100:500);
%set(ax2,'XTickLabel',{'0','','200','','400',''});
set(ax2,'XTickLabel',{'','100','','300','','500'});
set(ax2,'YTick',0:0.1:0.5);
set(ax2,'YTickLabel',[]);
Pk1 = plot(kActual,q,'k','Linewidth',1);
Pk2 = plot(kOpt,q,'k--','Linewidth',1);
%contour(x2,y2,z2,[1 1],'linewidth',0.5,'linestyle',':','color','k');
% shading interp
hold off
caxis([partition(1) partition(end)]); %fix the color limits to partition limits!
axis([0 500 0 0.5]);

%manual xlabels to bring them closer to axis
% make a vector of vertical positions after the offset:
%offset = repmat(ax2.YTick(1)-paper.xtickOff,1,numel(ax2.XTick));
% create new lables:
%text(ax2.XTick,offset,ax2.XTickLabel,'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex');
%ax2.XTickLabel=[]; %remove original ones

%set(ax2,'FontWeight','bold');
set(ax2,'FontSize',paper.fontSize);
%xl2 = xlabel('SE coupling $\lambda$','interpreter','latex','FontSize',paper.labelSize);
xl2 = xlabel(' no.~of qubits $k$','interpreter','latex','FontSize',paper.labelSize);
%xl2.Position(2) = -paper.xlabelOff;
ylabel([]);
tb = text(525,0.48,'(b)','interpreter','latex','FontSize',paper.labelSize+1);

%FUCK: text uses units relative to axis object, cannot set to colorbar!
Xfig2ax2=@(x)( (x-ax2.Position(1))/ax2.Position(3) );
Yfig2ax2=@(x)( (x-ax2.Position(2))/ax2.Position(4) );
tt = text(Xfig2ax2(hpos(1)+0.05)*500,Yfig2ax2(hpos(2)+hpos(4)+0.07)*0.5,'$\eta$','interpreter','latex','FontSize',paper.labelSize);
tt.Position


h=colorbar('eastoutside');
set(h,'Position',hpos);
set(h,'Ticks',linspace(0,0.8,10));
set(h,'TickLabels',{'0','','0.2','','0.4','','0.6','','0.8',''});
set(h,'FontSize',paper.fontSize);

figD.Units = 'centimeters';
figD.PaperUnits = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'TickLabelInterpreter','tex');
%figD.Position = [0 0 paper.width paper.height];
%figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[paper.width, paper.height]);
set(ax1, 'color', 'none');
set(ax2, 'color', 'none');
%set(h, 'color', 'none');
set(figD, 'color', 'none');
%print('fig1a.eps','-dpsc2');
print('fig5', '-dpdf', '-r600', '-loose')
close(figD);
