clf; 
% x = load('compareModels.dat');
x = [Sp,Sl,Sg]';
blue = [46 0 202]/255;
red = [227 0 115]/255;
plot(g,x(1,:),'color', blue,'linewidth',3);
hold on;

plot(g,x(3,:),'--','color',  blue,'linewidth',4);
plot(g,x(2,:),':', 'color', blue,'linewidth',4);
% 
% hold on;
% 
% semilogx(x(1,:),x(7,:)/0.01,'--','color',  red,'linewidth',4);
% semilogx(x(1,:),x(6,:)/0.01,':', 'color', red,'linewidth',4);
% 
% text(0.012,0.03,'Cooling','Rotation',90,'FontSize',16,'FontWeight','bold');
% text(0.012,-0.08,'Heating','Rotation',90,'FontSize',16,'FontWeight','bold');
axis([0.01 1 -0.02/100 0.04/100]);
grid off;
set(gca,'FontSize', 18);
set(gca,'YTick',(-0.02:0.02:0.05)/100);
% set(gca,'XTick',0:5:20);
xlabel('T_W-T_C');
ylabel('T_H-T_C');
print('test.eps','-dpsc2');