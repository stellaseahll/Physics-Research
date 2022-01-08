clf; clc; clear;
paper.width = 8;
paper.height = 5;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure(); 

% close(figure);

%% Plot Fisher Info for different K for 1 pass
% figD = figure();  

r =  1;
beta = logspace(-2,0,101);
Gamma = logspace(-2,0,101);
K = 0.1;
% K(end) = [];
p0 = 1;
% rho0 = diag([p0,1-p0]);
np = 1;
for i = 1:length(beta)
    i
    for j = 1:length(Gamma)
        FC(i,j)=getClassicalFishSteady(K,Gamma(j),beta(i),beta(i)*0.001,p0,np);
    end
end


contourf(log10(Gamma),log10(beta),FC);
xlabel('\bf Waiting time $\log_{10}(\Gamma)$ ','interpreter','latex','FontSize',8);
ylabel('\bf Inverse temperature $\log_{10}\beta$','interpreter','latex','FontSize',8);

figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig3b', '-dpdf', '-r600')
close(figD);