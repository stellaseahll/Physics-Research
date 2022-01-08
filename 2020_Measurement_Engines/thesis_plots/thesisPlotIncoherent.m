paper.width = 8;
paper.height = 2.5;
paper.fontSize = 12;
set(0,'defaulttextinterpreter','latex')
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

figD = figure(); 
c1 = [52 79 168]/255;
c2 = [193 24 0]/255;
c3 = [84 130 53]/255;
c4 = [207 138 68]/255;
c = [c1;c2;c3;c4];
% paper.axwidth = (paper.width-paper.leftgap-paper.N*paper.axgap)/paper.N;
% paper.axheight = paper.axwidth;
% paper.height = paper.bottomgap+paper.axheight+paper.axgap;

% paper.xtickOff = 0.65; %manual offset of xTickLabels from axis
% paper.xlabelOff = 0.09; %manual offset of xlabel from axis

% ax=subplot(1,2,1); 
a = 0.9;
fill([0 0 XEnd XEnd 0], [-0.1 1 1 -0.1 -0.1],[a a a],'linestyle','none');
hold on
plot([0 max(X)],[wth1 wth1],':','Color',c(1,:),'LineWidth',2);
plot([0 max(X)],[wth2 wth2],':','Color',c(2,:),'LineWidth',2);
plot([0 max(X)],[wth3 wth3],':','Color',c(3,:),'LineWidth',2);
plot([0 max(X)],[wth4 wth4],':','Color',c(4,:),'LineWidth',2);
plot(X,w4,'Color',c(4,:),'LineWidth',2);
plot(X,w3,'Color',c(3,:),'LineWidth',2);
plot(X,w2,'Color',c(2,:),'LineWidth',2);
plot(X,w1,'Color',c(1,:),'LineWidth',2);

xlabel('$x_0/x_{\mathrm th}$','interpreter','latex','FontSize',12);
ylabel('\bf Power $\mathcal{W}/\hbar\Omega\kappa_{\mathrm{h}}$','interpreter','latex','FontSize',12);
xlim([0 max(X)]);
%ylim([-0.02 0.9]);
ylim([0 1]);
set(gca,'XTick',[0,1,2,3],'XTickLabels',{'0','1','2','3'});
set(gca,'YTick',0:0.5:1.0);
set(gca,'FontWeight','bold');
% text(0.1,0.92,'(a)','interpreter','latex','FontSize',12);
% set(ax,'FontSize',paper.fontSize);

set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);

print('incoherentWork', '-dpdf', '-r600')
close(figD);

figD = figure();
fill([0 0 XEnd XEnd 0], [-0.1 1 1 -0.1 -0.1],[a a a],'linestyle','none');
hold on
plot([0 max(X)],[nth1 nth1],':','Color',c(1,:),'LineWidth',2);
plot([0 max(X)],[nth2 nth2],':','Color',c(2,:),'LineWidth',2);
plot([0 max(X)],[nth3 nth3],':','Color',c(3,:),'LineWidth',2);
plot([0 max(X)],[nth4 nth4],':','Color',c(4,:),'LineWidth',2);
plot(X,n4,'Color',c(4,:),'LineWidth',2);
plot(X,n3,'Color',c(3,:),'LineWidth',2);
plot(X,n2,'Color',c(2,:),'LineWidth',2);
plot(X,n1,'Color',c(1,:),'LineWidth',2);

% plot([0 max(X)],[1-wp/ws 1-wp/ws],'r','LineWidth',2.2); %Otto curve
x = logspace(-5,5,1000);
Tc = wp/2./acoth(x.^2);
plot(x0./x,1-Tc./Th,'k','LineWidth',2.2);%Carnot curve

xlabel('$x_0/x_{\rm th}$','interpreter','latex','FontSize',12);
ylabel('\bf Efficiency $\eta$','interpreter','latex','FontSize',12);
xlim([0 max(X)]);
ylim([0 1]);
set(gca,'XTick',[0,1,2,3],'XTickLabels',{'0','1','2','3'});
set(gca,'YTick',0:0.5:1.0);
set(gca,'FontWeight','bold');
% text(0.1,0.92,'(a)','interpreter','latex','FontSize',12);
% set(ax,'FontSize',paper.fontSize);

set(gca,'Units','centimeters');
pos = get(gca,'Position');
pos(3) = paper.width;
pos(4) = paper.height;
set(gca,'Position',pos);
%for saving
set(figD.Children,'FontSize',paper.fontSize);

print('incoherentEff', '-dpdf', '-r600')

close(figD);

