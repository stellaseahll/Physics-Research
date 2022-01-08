%% load data
panels = {'(a)','(b)'};

% load('8AncillaT2_g_e_p_tol6.mat');
% Fg = Fs{1}(1,:);
% Fe = Fs{2}(1,:);
% Fp = Fs{3}(1,:);

%here saved as arrays for a fixed gt=pi/100

% load('IncoherentDataNEW_wp1_ws100_g3.54_km0.1000_kh0.0010_kc0.1000.mat'); w1 = -Wq/ws/kh; n1 = real(-Wq)./real(Qh+(W-Wq)); n1(W>0) = NaN; wth1 = km*nh/kh/(2*nh+1+km/kh)*(1-g*g/ws/wp); nth1 = (1-g*g/ws/wp)/(1+(1+(2*nh+2)*kh/km)*g*g/ws/wp); 
% load('IncoherentDataNEW_wp1_ws100_g3.54_km0.0100_kh0.0010_kc0.1000.mat'); w2 = -Wq/ws/kh; n2 = real(-Wq)./real(Qh+(W-Wq)); n2(W>0) = NaN; wth2 = km*nh/kh/(2*nh+1+km/kh)*(1-g*g/ws/wp); nth2 = (1-g*g/ws/wp)/(1+(1+(2*nh+2)*kh/km)*g*g/ws/wp);
% load('IncoherentDataNEW_wp1_ws100_g3.54_km0.0010_kh0.0010_kc0.1000.mat'); w3 = -Wq/ws/kh; n3 = real(-Wq)./real(Qh+(W-Wq)); n3(W>0) = NaN; wth3 = km*nh/kh/(2*nh+1+km/kh)*(1-g*g/ws/wp); nth3 = (1-g*g/ws/wp)/(1+(1+(2*nh+2)*kh/km)*g*g/ws/wp);
% load('IncoherentDataNEW_wp1_ws100_g3.54_km0.0001_kh0.0010_kc0.1000.mat'); w4 = -Wq/ws/kh; n4 = real(-Wq)./real(Qh+(W-Wq)); n4(W>0) = NaN; wth4 = km*nh/kh/(2*nh+1+km/kh)*(1-g*g/ws/wp); nth4 = (1-g*g/ws/wp)/(1+(1+(2*nh+2)*kh/km)*g*g/ws/wp);

%For consistency check of old and new computations
load('IncoherentDataNEW_wp1_ws100_g3.54_km0.1000_kh0.0010_kc0.1000.mat'); 
w1b = -real(W)/ws/kh; n1b = real(-W)./real(Qh); n1b(W>0) = NaN; 
load('IncoherentDataNEW_wp1_ws100_g3.54_km0.0100_kh0.0010_kc0.1000.mat'); 
w2b = -real(W)/ws/kh; n2b = real(-W)./real(Qh); n2b(W>0) = NaN; 
load('IncoherentDataNEW_wp1_ws100_g3.54_km0.0010_kh0.0010_kc0.1000.mat'); 
w3b = -real(W)/ws/kh; n3b = real(-W)./real(Qh); n3b(W>0) = NaN; 
load('IncoherentDataNEW_wp1_ws100_g3.54_km0.0001_kh0.0010_kc0.1000.mat'); 
w4b = -real(W)/ws/kh; n4b = real(-W)./real(Qh); n4b(W>0) = NaN; 

%New computations to be used
load('IncoherentDataNEW2_wp1_ws100_g3.54_km0.1000_kh0.0010_kc0.1000.mat'); 
w1 = -real(W)/ws/kh; n1 = real(-W)./real(Qh); n1(W>0) = NaN; 
wth1 = km*nh/kh/(2*nh+1+km/kh)*(1-2*wp*x0^2/ws); 
nth1 = (1-2*wp*x0^2/ws)/(1+(1+(2*nh+2)*kh/km)*2*wp*x0^2/ws); 
qc1 = -real(Qc)./Tc/log(2)/km;
s1 = -real(Qc./Tc+Qh./Th);
load('IncoherentDataNEW2_wp1_ws100_g3.54_km0.0100_kh0.0010_kc0.1000.mat'); 
w2 = -real(W)/ws/kh; n2 = real(-W)./real(Qh); n2(W>0) = NaN; 
wth2 = km*nh/kh/(2*nh+1+km/kh)*(1-2*wp*x0^2/ws); 
nth2 = (1-2*wp*x0^2/ws)/(1+(1+(2*nh+2)*kh/km)*2*wp*x0^2/ws);
qc2 = -real(Qc)./Tc/log(2)/km;
s2 = -real(Qc./Tc+Qh./Th);
load('IncoherentDataNEW2_wp1_ws100_g3.54_km0.0010_kh0.0010_kc0.1000.mat'); 
w3 = -real(W)/ws/kh; n3 = real(-W)./real(Qh); n3(W>0) = NaN; 
wth3 = km*nh/kh/(2*nh+1+km/kh)*(1-2*wp*x0^2/ws); 
nth3 = (1-2*wp*x0^2/ws)/(1+(1+(2*nh+2)*kh/km)*2*wp*x0^2/ws);
qc3 = -real(Qc)./Tc/log(2)/km;
s3 = -real(Qc./Tc+Qh./Th);
load('IncoherentDataNEW2_wp1_ws100_g3.54_km0.0001_kh0.0010_kc0.1000.mat'); 
w4 = -real(W)/ws/kh; n4 = real(-W)./real(Qh); n4(W>0) = NaN; 
wth4 = km*nh/kh/(2*nh+1+km/kh)*(1-2*wp*x0^2/ws); 
nth4 = (1-2*wp*x0^2/ws)/(1+(1+(2*nh+2)*kh/km)*2*wp*x0^2/ws);
qc4 = -real(Qc)./Tc/log(2)/km;
s4 = -real(Qc./Tc+Qh./Th);

%Tc = wp./log(1./nc+1);
%xth = sqrt(coth(wp/2./Tc));
% X = (g/wp)./xth;
%x0 = g/wp/sqrt(2);
X = x0./xth;
TcEnd = wp./log(1+1);
xthEnd = sqrt(coth(wp/2./TcEnd));
XEnd = x0./xthEnd;

%% plot wide fig (small g) 

paper.leftgap = 1.3;
paper.bottomgap = 0.9;
paper.axgap = 0.5;
paper.width = 14;
paper.height = 6;
paper.N = 2;
set(0,'defaulttextinterpreter','latex');
figD = figure(666);
c = colormap(lines);
% paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
% paper.axheight = paper.axwidth;
% paper.height = paper.bottomgap+paper.axheight+paper.axgap;

paper.fontSize = 12;
paper.labelSize = 12;
% paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
% paper.xlabelOff = 0.09; %manual offset of xlabel from axis

ax=subplot(1,2,1); 
a = 0.9;
fill([0 0 XEnd XEnd 0], [-0.1 1 1 -0.1 -0.1],[a a a],'linestyle','none');
hold on
plot(X,w1,'Color',c(1,:),'LineWidth',1);
plot(X,w2,'Color',c(2,:),'LineWidth',1);
plot(X,w3,'Color',c(4,:),'LineWidth',1);
plot(X,w4,'Color',c(5,:),'LineWidth',1);
plot([0 max(X)],[wth1 wth1],':','Color',c(1,:),'LineWidth',1);
plot([0 max(X)],[wth2 wth2],':','Color',c(2,:),'LineWidth',1);
plot([0 max(X)],[wth3 wth3],':','Color',c(4,:),'LineWidth',1);
plot([0 max(X)],[wth4 wth4],':','Color',c(5,:),'LineWidth',1);
xlabel('$x_0/x_{\rm th}$','interpreter','latex','FontSize',paper.labelSize);
ylabel('Power $\dot{W}/\hbar\Omega\kappa_h$','interpreter','latex','FontSize',paper.labelSize);
xlim([0 max(X)]);
%ylim([-0.02 0.9]);
ylim([0 1]);
set(ax,'XTick',[1,2,3],'XTickLabels',{'1','2','3'});
set(ax,'YTick',0:0.2:1.0);
text(0.1,0.92,'(a)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);

ax = subplot(1,2,2);
fill([0 0 XEnd XEnd 0], [-0.1 1 1 -0.1 -0.1],[a a a],'linestyle','none');
hold on
plot(X,n1,'Color',c(1,:),'LineWidth',1);
plot(X,n2,'Color',c(2,:),'LineWidth',1);
plot(X,n3,'Color',c(4,:),'LineWidth',1);
plot(X,n4,'Color',c(5,:),'LineWidth',1);
plot([0 max(X)],[nth1 nth1],':','Color',c(1,:),'LineWidth',1);
plot([0 max(X)],[nth2 nth2],':','Color',c(2,:),'LineWidth',1);
plot([0 max(X)],[nth3 nth3],':','Color',c(4,:),'LineWidth',1);
plot([0 max(X)],[nth4 nth4],':','Color',c(5,:),'LineWidth',1);
% plot([0 max(X)],[1-wp/ws 1-wp/ws],'r','LineWidth',1.2); %Otto curve
x = logspace(-5,5,1000);
Tc = wp/2./acoth(x.^2);
plot(x0./x,1-Tc./Th,'k','LineWidth',1.2);%Carnot curve

xlabel('$x_0/x_{\rm th}$','interpreter','latex','FontSize',paper.labelSize);
ylabel('Efficiency $\eta$','interpreter','latex','FontSize',paper.labelSize);
xlim([0 max(X)]);
ylim([0 1]);
set(ax,'XTick',1:3);
set(ax,'YTick',0:0.2:1);
text(0.1,0.92,'(b)','interpreter','latex','FontSize',paper.labelSize);
set(ax,'FontSize',paper.fontSize);


% 
figD.Units = 'centimeters';
figD.Position(3) = paper.width;
figD.Position(4) = paper.height;
set(figD.Children,'TickLabelInterpreter','latex');
figD.PaperPositionMode = 'auto';
posD=get(figD,'Position');
set(figD,'PaperSize',[posD(3), posD(4)]);
set(figD, 'color', 'none');
% %print('fig2.eps','-dpsc2');
print('fig.pdf', '-dpdf', '-r600')
% %fix_lines('figF2.eps','figF2b.eps');
%system('/Library/TeX/texbin/epstopdf figF2b.eps')
%export_fig figF2b.eps
close(figD);

