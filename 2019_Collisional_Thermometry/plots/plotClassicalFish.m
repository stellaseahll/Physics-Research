clf; clc; clear;
paper.width = 8;
paper.height = 3;
paper.fontSize = 8;
set(0,'defaulttextinterpreter','latex');
c = colormap(lines);
figD = figure(); 

% close(figure);

%% Plot Fisher Info for different K for 1 pass
% figD = figure();  

r =  1;
Gamma = 1;
np = 5:5:20;
beta = 0.01:0.01:1;
K = 0;
for i = 1:length(beta)
    F1=getClassicalFish(K,Gamma,beta(i),beta(i)*0.001,r,10,1);
    for j = 1:length(np)
        F(i,j)=getClassicalFish(K,Gamma,beta(i),beta(i)*0.001,r,np(j));
        Fnorm(i,j) =  F(i,j)/np(j)/F1;
    end
    
end

plot(beta,Fnorm','LineWidth',1);
set(gca,'FontWeight','bold');
xlabel('\bf inverse temperature $\beta$ ','interpreter','latex','FontSize',8);
ylabel('\bf $F_N/NF_1$','interpreter','latex','FontSize',8);

%for saving
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig1', '-dpdf', '-r600')
close(figD);
% 
