
set(0,'defaulttextinterpreter','latex');
%figX.Units = 'centimeters';
%figX.PaperUnits = 'centimeters';

figX = figure(667); 

% load('fig3_unbounded.mat');
tend = 500;
myred = [198 86 106]/255;
% c(1,:) = [52 79 168]/255;
% c(2,:) = [193 24 0]/255;
% c(3,:) = [84 130 53]/255;
% c(4,:) = [207 138 68]/255;
% c(5,:) = [0 0 0]/255;
% runInfLadderData;

load('fig3NEW_q0.25ground.mat');
EI = Einc;
EC = Ecoh;
ErgI = Erginc;
ErgC = Ergcoh;
ErgD = Ergdeph;
load('fig3NEW_q0.49ground.mat');
EI(2,:) = Einc;
EC(2,:) = Ecoh;
ErgI(2,:) = Erginc;
ErgC(2,:) = Ergcoh;
ErgD(2,:) = Ergdeph;


%% Contour plot

%fileInName = {'fig3NEW_q0.25ground.mat','fig3NEW_q0.49ground.mat'};
fileOutName = {'fig3NEW_q025ground','fig3NEW_q049ground'};
labelName = {'(a)','(b)'};

tend = 2500;
paper.width = 8;
paper.height = 3.5;
paper.bottomgap = 0.7;
paper.fontSize = 7;
paper.labelSize = 8;
paper.bottomgap = 0.75;

axpos = [1;1]*[0.84/paper.width, paper.bottomgap/paper.height, 3.45/paper.width, 2.6/paper.height];
axpos(2,1) = 4.45/paper.width;
xpos = [292 392];

for r = 1:2
    %load(fileInName{r});
    % c(1,:) = [52 79 168]/255;
    % c(2,:) = [193 24 0]/255;
    % c(3,:) = [84 130 53]/255;
    % c(4,:) = [207 138 68]/255;
    % c(5,:) = [0 0 0]/255;
    % runInfLadderData;

    a(r) = subplot('position',axpos(r,:));
    plot(0:tend,EI(r,1:tend+1),'k','linewidth',1);
    hold on;
    plot(0:tend,EC(r,1:tend+1),'color',myred,'linewidth',1);
    plot(0:tend,ErgI(r,1:tend+1),'k--','linewidth',1);
    plot(0:tend,ErgC(r,1:tend+1),'color',myred,'linestyle','--','linewidth',1);
    plot(0:tend,ErgD(r,1:tend+1),'color',myred,'linestyle',':','linewidth',1);
    plot(xpos(r)*[1 1],[0 200],'color',[0.7 0.7 0.7],'linestyle',':','linewidth',1);

    hold off
    set(a(r),'YTick',0:100:N);
    set(a(r),'XTick',0:1000:tend);
    set(a(r),'FontSize',paper.fontSize);
    if (r>1)
        set(a(r),'YTickLabel',[]);
    else
        ylabel('energy in units $E$','interpreter','latex','FontSize',paper.fontSize+1);
    end
    % set(gca,'YTick',0:20:60);
    axis([0 tend 0 N]);
    % R = [ (sE(1:tend+1)-sE(1)) (sErg(1:tend+1)-sErg(1)) (tE(1:tend+1)-tE(1)) (tErg(1:tend+1)-tErg(1))  (tDeErg(1:tend+1)-tDeErg(1))];

    text(2150,0.9*N,labelName{r},'FontSize',paper.fontSize+2);
    xlabel('no. of qubits $k$','interpreter','latex','FontSize',paper.fontSize+1);
    % print(fileOutName{i}, '-dpdf', '-r600');
    % close(figD);
end

figX.Units = 'centimeters';
figX.PaperUnits = 'centimeters';
figX.Position(3) = paper.width;
figX.Position(4) = paper.height;
%set(figX.Children,'FontSize',paper.fontSize);
set(figX, 'color', 'none');
set(figX,'PaperSize',[paper.width, paper.height]);
figX.PaperPositionMode = 'auto';
print('fig3NEW', '-dpdf', '-r600', '-loose');
%export_fig fig3.pdf -painters
%my_fix_lines('fig3NEW.pdf','fig3.pdf');
close(figX);