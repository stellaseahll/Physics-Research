%% load data
panels = {'(a)','(b)','(c)','(d)'};

% load('8AncillaT2_g_e_p_tol6.mat');
% Fg = Fs{1}(1,:);
% Fe = Fs{2}(1,:);
% Fp = Fs{3}(1,:);

T = 2; %this is kB*T/(hbar*omega)
nbar = 1/(exp(1/T)-1);
p1th = nbar/(2*nbar+1);
Fth = 1 / ( (exp(1/T)/(exp(1/T)-1)^2)^2/(p1th-p1th^2) );

%here saved as arrays for a fixed gt=pi/100
load('12AncillaT2_State1_gt0.01.mat'); Fg = F;
load('12AncillaT2_State2_gt0.01.mat'); Fe = F;
load('12AncillaT2_State3_gt0.01.mat'); Fp = F;

%here saved as arrays for a fixed gt=pi/4
load('12AncillaT2_State1_gt0.25.mat'); Fg2 = F;
load('12AncillaT2_State2_gt0.25.mat'); Fe2 = F;
load('12AncillaT2_State3_gt0.25.mat'); Fp2 = F;

np = (1:12);
%% plot wide fig (small g) 

paper.leftgap = 1.3;
paper.bottomgap = 0.9;
paper.axgap = 0.3;
paper.width = 16;
paper.height = 5;
paper.N = 3;

paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
paper.axheight = paper.axwidth;
paper.height = paper.bottomgap+paper.axheight+paper.axgap;

paper.fontSize = 10;
paper.labelSize = 10;
paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
paper.xlabelOff = 0.09; %manual offset of xlabel from axis

set(0,'defaulttextinterpreter','latex');
figD = figure(666);
c = colormap(lines);
for n=1:paper.N

    ax{n} = subplot('position', ...
        [ (paper.leftgap+(n-1)*(paper.axgap+paper.axwidth))/paper.width, ...
        paper.bottomgap/paper.height, paper.axwidth/paper.width, paper.axheight/paper.height]);
% old cell scheme of stored F-values
%     loglog(ax{n},np,Fg{n}(1)*np/Fth,':','Color',c(1,:),'LineWidth',0.5)
%     hold on
%     loglog(ax{n},np,Fe{n}(1)*np/Fth,':','Color',c(2,:),'LineWidth',0.5)
%     loglog(ax{n},np,Fp{n}(1)*np/Fth,':','Color',c(4,:),'LineWidth',0.5)
%     loglog(ax{n},np,Fg{n}/Fth,'.','Color',c(1,:),'MarkerSize',10)
%     loglog(ax{n},np,Fe{n}/Fth,'o','Color',c(2,:),'MarkerSize',5)
%     loglog(ax{n},np,Fp{n}/Fth,'+','Color',c(4,:),'MarkerSize',6)
    loglog(ax{n},np,Fg(n,1)*np/Fth,':','Color',c(1,:),'LineWidth',0.5)
    hold on
    loglog(ax{n},np,Fe(n,1)*np/Fth,':','Color',c(2,:),'LineWidth',0.5)
    loglog(ax{n},np,Fp(n,1)*np/Fth,':','Color',c(4,:),'LineWidth',0.5)
    loglog(ax{n},np,Fg(n,:)/Fth,'.','Color',c(1,:),'MarkerSize',10)
    loglog(ax{n},np,Fe(n,:)/Fth,'o','Color',c(2,:),'MarkerSize',5)
    loglog(ax{n},np,Fp(n,:)/Fth,'+','Color',c(4,:),'MarkerSize',6)
    xlim([0.9 14]);
    ylim([8e-5 6e-2]);
    hold off
    set(ax{n},'XTick',[1,5,10],'XTickLabels',{'1','5','10'});
    set(ax{n},'FontSize',paper.fontSize);
    ax{n}.XAxis.MinorTickValues=1:12;
    %manual xlabels to bring them closer to axis
    % make a vector of vertical positions after the offset:
    offset = repmat(5.6e-5,1,numel(ax{n}.XTick));
    % create new lables:
    text(ax{n}.XTick,offset,ax{n}.XTickLabel,'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex');
    xl{n} = xlabel('$N$','interpreter','latex','FontSize',paper.labelSize);
    xl{n}.Position(2) = 4.5e-5;
    ax{n}.XTickLabel=[]; %remove original ones
    
    %panel label
    text(1.02,3.1e-2,panels{n},'interpreter','latex','FontSize',paper.labelSize)

    if n==1 %only leftmost ylabel
        yl = ylabel('$\mathcal{F}_N/\mathcal{F}_{\rm th}$','interpreter','latex','FontSize',paper.labelSize);
        yl.Position(1) = 0.57;
    else
        ax{n}.YTickLabel=[]; %no ticks
    end
    
end

figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'TickLabelInterpreter','latex');
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
%print('fig2.eps','-dpsc2');
print('fig2.pdf', '-dpdf', '-r600')
%fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);

%% plot wide fig (large g)

paper.leftgap = 1.3;
paper.bottomgap = 0.9;
paper.axgap = 0.3;
paper.width = 16;
paper.height = 5;
paper.N = 3;

paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
paper.axheight = paper.axwidth;
paper.height = paper.bottomgap+paper.axheight+paper.axgap;

paper.fontSize = 10;
paper.labelSize = 10;
paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
paper.xlabelOff = 0.09; %manual offset of xlabel from axis

set(0,'defaulttextinterpreter','latex');
figD = figure(666);
c = colormap(lines);

for n=1:paper.N

    ax{n} = subplot('position', ...
        [ (paper.leftgap+(n-1)*(paper.axgap+paper.axwidth))/paper.width, ...
        paper.bottomgap/paper.height, paper.axwidth/paper.width, paper.axheight/paper.height]);
% old cell scheme of stored F-values
%     loglog(ax{n},np,Fg{n}(1)*np/Fth,':','Color',c(1,:),'LineWidth',0.5)
%     hold on
%     loglog(ax{n},np,Fe{n}(1)*np/Fth,':','Color',c(2,:),'LineWidth',0.5)
%     loglog(ax{n},np,Fp{n}(1)*np/Fth,':','Color',c(4,:),'LineWidth',0.5)
%     loglog(ax{n},np,Fg{n}/Fth,'.','Color',c(1,:),'MarkerSize',10)
%     loglog(ax{n},np,Fe{n}/Fth,'o','Color',c(2,:),'MarkerSize',5)
%     loglog(ax{n},np,Fp{n}/Fth,'+','Color',c(4,:),'MarkerSize',6)
    loglog(ax{n},np,Fg2(n,1)*np/Fth,':','Color',c(1,:),'LineWidth',0.5)
    hold on
    loglog(ax{n},np,Fe2(n,1)*np/Fth,':','Color',c(2,:),'LineWidth',0.5)
    loglog(ax{n},np,Fp2(n,1)*np/Fth,':','Color',c(4,:),'LineWidth',0.5)
    loglog(ax{n},np,Fg2(n,:)/Fth,'.','Color',c(1,:),'MarkerSize',10)
    loglog(ax{n},np,Fe2(n,:)/Fth,'o','Color',c(2,:),'MarkerSize',5)
    loglog(ax{n},np,Fp2(n,:)/Fth,'+','Color',c(4,:),'MarkerSize',6)
    xlim([0.9 14]);
    ylim([1.5e-3 5e1]);
    hold off
    set(ax{n},'XTick',[1,5,10],'XTickLabels',{'1','5','10'},'YTick',[1e-2,1e-1,1e0,1e1]);
    set(ax{n},'FontSize',paper.fontSize);
    ax{n}.XAxis.MinorTickValues=1:12;
    %manual xlabels to bring them closer to axis
    % make a vector of vertical positions after the offset:
    offset = repmat(8e-4,1,numel(ax{n}.XTick));
    % create new lables:
    text(ax{n}.XTick,offset,ax{n}.XTickLabel,'HorizontalAlign','center','FontSize',paper.fontSize,'interpreter','latex');
    xl{n} = xlabel('$N$','interpreter','latex','FontSize',paper.labelSize);
    xl{n}.Position(2) = 5.8e-4;
    ax{n}.XTickLabel=[]; %remove original ones
    
    %panel label
    text(1.02,21,panels{n},'interpreter','latex','FontSize',paper.labelSize)

    if n==1 %only leftmost ylabel
        yl = ylabel('$\mathcal{F}_N/\mathcal{F}_{\rm th}$','interpreter','latex','FontSize',paper.labelSize);
        yl.Position(1) = 0.55;
    else
        ax{n}.YTickLabel=[]; %no ticks
    end
    
end

figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'TickLabelInterpreter','latex');
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
%print('fig2_halfswap.eps','-dpsc2');
print('fig2_halfswap.pdf', '-dpdf', '-r600')
%fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);