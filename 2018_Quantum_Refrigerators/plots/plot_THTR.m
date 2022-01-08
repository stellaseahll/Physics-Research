clear;
set(0,'defaulttextinterpreter','latex');
paper.width = 8;
paper.height = 2.5;
paper.fontSize = 12;

THdata = load('THdata0.01.dat');
TW = THdata(1,:)+1;
THcool = 5./(1/1+4./TW);
ymax = THcool(end)+0.1;
hot = [TW TW(end) TW(1) TW(1); THcool ymax ymax 1]-1;
cold = [TW TW(end) TW(1); THcool 0 0]-1;
c1 = [0 82 177]/255;
c2 = [192 0 0]/255;
c3 = [84 130 53]/255;
c4 = [207 138 68]/255;
red = c2;
blue = c1;
figD = figure;
% fill(hot(1,:),hot(2,:),red,'edgecolor','none');
fill(hot(1,:),hot(2,:),[220 220 220]/255,'edgecolor','none');hold on;
% fill(cold(1,:),cold(2,:),blue,'edgecolor','none');
% fill([0 (TW-1) 20 0],[0 THcool-1 0 0],[220 220 220]/255);hold on;
plot(TW-1,THcool-1,'k','linewidth',2);
% plot(THdata(1,:),THdata(2:8,:),'k','linewidth',2);
col = linspace(0.8,0.3,4);
R=1:2:8;
for k = 1:4
    i = R(k);
    plot(THdata(1,:),THdata(10-i,:),'color',[1 1 1]*col(k),'linewidth',1);
%     plot(THdata(1,1:i*10-2),THdata(10-i,1:i*10-2),'k','linewidth',2);
%     plot(THdata(1,i*10+7:end),THdata(10-i,i*10+7:end),'k','linewidth',2);
%     text(THdata(1,i*10),THdata(10-i,i*10),sprintf('%.1f',0.9-i*0.1),'FontSize',12,'FontWeight','bold');
end
% 
% for i = 1
%     plot(THdata(1,1:i*10-2),THdata(10-i,1:i*10-2),'k','linewidth',2);
%     plot(THdata(1,i*10+4:end),THdata(10-i,i*10+4:end),'k','linewidth',2);
%     text(THdata(1,i*10),THdata(10-i,i*10),sprintf('%.1f',0.9-i*0.1),'FontSize',12,'FontWeight','bold');
% end
% text(THdata(1,end)+0.5,THdata(2,end),'0.1','FontSize',20);
% text(THdata(1,end)+0.5,THdata(3,end),'0.2','FontSize',20);
% text(THdata(1,end)+0.5,THdata(4,end),'0.3','FontSize',20);
% text(THdata(1,end)+0.5,THdata(5,end),'0.4','FontSize',20);
% text(THdata(1,50)+2,THdata(6,50),'0.5','FontSize',20);
% text(THdata(1,40)-1,THdata(7,50)+0.54,'0.6','FontSize',20);
% text(THdata(1,30)-2,THdata(8,20)+0.1,'0.7','FontSize',20);
axis([0 20 0 ymax-1]);
% set(gca,'FontSize', 18);
set(gca,'YTick',[0:3]);
set(gca,'XTick',0:5:20);
xlabel('$T_{\mathrm {w}}-T_{\mathrm {c}}$','interpreter','latex','FontSize',12);
ylabel('$T_{\mathrm {h}}-T_{\mathrm {c}}$','interpreter','latex','FontSize',12);
set(gca,'FontWeight','bold');
% set(gca, 'XTickLabel', []);
set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
set(figD.Children,'FontSize',paper.fontSize);
print('coolingwindow', '-dpdf', '-r600');
close(figD);
% xlabel('T_W-T_C');
% ylabel('T_H-T_C');
% print('test.eps','-dpsc2');