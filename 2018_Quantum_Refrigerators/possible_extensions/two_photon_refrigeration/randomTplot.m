%%Script for producing scatter plot of maximum QC values and g and chi
%%given random temperatures

clear;

EC = 1;
EW = 1;
EH = EC+EW;

load(['maxQCrandomT_Ec1_Ew',num2str(EW),'_N5000_sqp.mat']);
snam = ['randTempEw',num2str(EW),'.pdf'];
lab = ['$\omega_w/\omega_c = ',num2str(EW),'$'];


kappa = chi/6 .* ( 2*EC./(exp(EC./Tc)-1) + 2*EW./(exp(EW./Tw)-1) + 2*EH./(exp(EH./Th)-1) +  EH*2);

%select subset of sample (QC>0 and chi<0.1)
w = (QC>=0)&(chi<0.1);
QC = QC(w);
chi = chi(w);
g = g(w);
kappa = kappa(w);

%sorting from small QC to large
[~,I] = sort(QC);
QC = QC(I);
chi = chi(I);
g = g(I);
kappa = kappa(I);

min(g)

[min(kappa),max(kappa)]


%% PLOT

paper.width = 5;
paper.height = 4;
paper.fontSize = 7;
set(0,'defaulttextinterpreter','latex')

blue = [46 0 202]/255;
red = [227 0 115]/255;

%color vector
CC = red.'*ones(1,length(QC)) + (blue-red).' * (QC-min(QC))/(max(QC)-min(QC));
%size vector
s = 15*ones(1,length(QC)) - 10*(QC-min(QC))/(max(QC)-min(QC));

figD = figure();
loglog([0.001, 1],[0.001,1],'Color',0.7*[1 1 1],'LineWidth',1)
xlim(gca,[0.002 1])
ylim(gca, [0.002 1])
set(gca,'YTick',[1e-2,1e-1,1]);
set(gca,'YTickLabel',{'0.01','0.1','1'});
set(gca,'XTick',[1e-3,1e-2,1e-1,1]);
set(gca,'XTickLabel',{'','0.01','0.1','1'});
set(gca,'Layer','top');
%set(gca,'XGrid','on')
%grid on
%grid minor
hold on
scatter(g,kappa,s.',CC.','.')
%grid on
%set(gca,'xscale','log')
%set(gca,'yscale','log')

hold off
set(figD.Children,'FontSize',paper.fontSize,'FontWeight','bold')

xlabel('Coupling strength $g/\omega_c$','interpreter','latex','FontSize',10)
ylabel('Therm.~rate $\kappa_{\rm eff}/\omega_c$','interpreter','latex','FontSize',10)

text(0.05,0.0040,lab,'FontSize',10,'interpreter','latex');


%axis([0.01 1 -0.11 0.19]);
%set(gca,'FontSize', 18);
%xlabel('T_W-T_C');
%ylabel('T_H-T_C');



%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
figD.PaperPositionMode = 'auto';
figD.Renderer='Painters';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)])
set(figD, 'color', 'none');
print(snam, '-dpdf', '-r600')
close(figD)