load('fig2_unbounded.mat');

set(0,'defaulttextinterpreter','latex');
paper.width = 8;
paper.height = 3.5;
paper.bottomgap = 0.7;
paper.fontSize = 7;
paper.labelSize = 8;
paper.bottomgap = 0.75;

myred = [198 86 106]/255;

axpos = [1;1]*[0.84/paper.width, paper.bottomgap/paper.height, 3.25/paper.width, 2.73/paper.height];
axpos(2,1) = 4.5/paper.width;
lw = 0.8;
lsv = '--'
lwv = 0.5;


int = 100;
kend = 501;
load('fig2_unboundedInf.mat');
theta = pi/2;
q = 0.25;
c = sqrt(q*(1-q));
ptheta = sin(theta/2)^2;

nR = 50-ptheta*(2*q-1)*(0:tend)+2*c*sqrt(ptheta*(1-ptheta))*(0:tend);
nL = 50-ptheta*(2*q-1)*(0:tend)-2*c*sqrt(ptheta*(1-ptheta))*(0:tend);
nM = 50-ptheta*(2*q-1)*(0:tend);
k = 1:int:kend;
distHeight = int*0.9;
n0 = find(sDist(1,:)==1);
xidx = n0-50:n0+150;
figD = figure(); 
subplot('position',axpos(1,:));
for i = 1:length(k)
    distx = k(i)-0.2+sDist(k(i),xidx)/max(sDist(k(i),xidx))*distHeight;
    disty = (k(i))+tDist(k(i),xidx)/max(tDist(k(i),xidx))*distHeight;
    plot(xidx-n0+50,distx,'color','k','linewidth',lw); hold on;
    plot(xidx-n0+50,disty,'color',myred,'linewidth',lw); 
    h=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50],[distx k(i)-0.2 k(i)-0.2 distx(1)],'k'); set(h,'facealpha',0.5);set(h,'linestyle','none');
    b=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50],[disty k(i)-0.2 k(i)-0.2 disty(1)],myred); set(b,'facealpha',0.5);set(b,'linestyle','none')
    plot([nR(k(i)) nR(k(i))],[k(i) k(i)+distHeight*1.2],'color',myred,'linewidth',lwv,'linestyle',lsv); 
    plot([nL(k(i)) nL(k(i))],[k(i) k(i)+distHeight*1.2],'color',myred,'linewidth',lwv,'linestyle',lsv);
    plot([nM(k(i)) nM(k(i))],[k(i) k(i)+distHeight*1.2],'color','k','linewidth',lwv,'linestyle',lsv);
end
distx = k(end)+int-0.2+sDistInf/max(sDistInf)*distHeight;
disty = k(end)+int-0.2+tDistInf/max(tDistInf)*distHeight;
plot(xidx-n0+50,distx,'color','k','linewidth',lw);
plot(xidx-n0+50,disty,'color',myred,'linewidth',lw);
h=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50],[distx k(end)+int-0.2 k(end)+int-0.2 distx(1)],'k'); set(h,'facealpha',0.5);set(h,'linestyle','none');
b=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50],[disty k(end)+int-0.2 k(end)+int-0.2 disty(1)],myred); set(b,'facealpha',0.5);set(b,'linestyle','none')
axis([0 200 0 kend+int-1]);
text(170,48,'(a)','FontSize',paper.fontSize+2);
% axis([0 25 0 0.14]);
set(gca,'YTick',0:int:kend);

set(gca,'XTick',0:50:249);
set(gca,'FontSize',paper.fontSize);
ylabel('no. of qubits $k$','interpreter','latex','FontSize',paper.fontSize+1);
xlabel('energy level $n$','interpreter','latex','FontSize',paper.fontSize+1);

load('fig2_bounded.mat');
load('fig2_boundedInf.mat');

subplot('position',axpos(2,:));

k = 1:int:kend;
distHeight = int*0.9;
n0 = find(sDist(1,:)==1);
xidx = n0:n0+200;
for i = 1:length(k)
    distx = k(i)-0.2+sDist(k(i),xidx)/max(sDist(k(i),xidx))*distHeight;
    disty = (k(i))+tDist(k(i),xidx)/max(tDist(k(i),xidx))*distHeight;
    plot(xidx-n0,distx,'color','k','linewidth',lw); hold on;
    plot(xidx-n0,disty,'color',myred,'linewidth',lw); 
    h=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50]-50,[distx k(i)-0.2 k(i)-0.2 distx(1)],'k'); set(h,'facealpha',0.5);set(h,'linestyle','none');
    b=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50]-50,[disty k(i)-0.2 k(i)-0.2 disty(1)],myred); set(b,'facealpha',0.5);set(b,'linestyle','none')
     plot([nR(k(i)) nR(k(i))]-50,[k(i) k(i)+distHeight*1.2],'color',myred,'linewidth',lwv,'linestyle',lsv); 
    plot([nL(k(i)) nL(k(i))]-50,[k(i) k(i)+distHeight*1.2],'color',myred,'linewidth',lwv,'linestyle',lsv);
    plot([nM(k(i)) nM(k(i))]-50,[k(i) k(i)+distHeight*1.2],'color','k','linewidth',lwv,'linestyle',lsv);
end
distx = k(end)+int-0.2+sDistInf/max(sDistInf)*distHeight;
disty = k(end)+int-0.2+tDistInf/max(tDistInf)*distHeight;
plot(xidx-n0,distx,'color','k','linewidth',lw);
plot(xidx-n0,disty,'color',myred,'linewidth',lw);
h=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50]-50,[distx k(end)+int-0.2 k(end)+int-0.2 distx(1)],'k'); set(h,'facealpha',0.5);set(h,'linestyle','none');
b=fill([xidx-n0+50 xidx(end)-n0+50 xidx(1)-n0+50 xidx(1)-n0+50]-50,[disty k(end)+int-0.2 k(end)+int-0.2 disty(1)],myred); set(b,'facealpha',0.5);set(b,'linestyle','none')

axis([0 200 0 kend+int-1]);
text(170,48,'(b)','FontSize',paper.fontSize+2);% axis([0 25 0 0.14]);
set(gca,'YTick',[]);
set(gca,'XTick',-50:50:249);
set(gca,'FontSize',paper.fontSize);
% ylabel('no. of qubits $k$','interpreter','latex','FontSize',12);
xlabel('energy level $n$','interpreter','latex','FontSize',paper.fontSize+1);



figD.Units = 'centimeters';
figD.PaperUnits = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
%set(figD.Children,'FontSize',paper.fontSize);
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
print('fig2NEW', '-dpdf', '-r600', '-loose');
close(figD);