%compare coherent VS incoherent as function of interrogation rates

clear;
wp = 1;
ws = 100;
% nc = 1;
kc = 0.1;
kh = 0.001;
%g = 2.5;
G = 0.1;
nh = 1.0;
Th = ws/log(1/nh+1);
dissType = 'g';

nc = 1.0;
Tc = wp./log(1./nc+1);
xth = sqrt(coth(wp/2/Tc));

x0 = 2.5;
g = sqrt(2)*wp*x0;

delta = 2*x0^2*wp;

NI = 201;
NC = 151;

%km = (1:N)/N * 0.1*wp;
km = wp*logspace(-4,0,NI);
G = wp*logspace(-3,0,NC);

%% coherent run 

dimpC = 45;
% 
% f=1;
% for i = 1:N
%     tic
%     fprintf('Coherent flat, step %i -- ',i);
%     M = modelCoherentEngine(dimpC,ws,wp,g,G(i),kh,kc,Th,Tc,f,dissType,delta);
%     M.findSS;
%     [WC1(i), QhC1(i), QcC1(i)] = M.getThermoProp;
%     toc
%     clear M;
%     %filename = sprintf('CoherentDataNEW3_G_Global_wp%.0f_ws%.0f_g%.2f_kh%.4f_kc%.4f.mat',wp,ws,g,kh,kc);
%     %save(filename);
% end

f = @(x) (1-sign(x))/2;
for i = 1:NC
    tic
    fprintf('Coherent step, step %i -- ',i);
    M = modelCoherentEngine(dimpC,ws,wp,g,G(i),kh,kc,Th,Tc,f,dissType,delta);
    M.findSS;
    [WC(i), QhC(i), QcC(i)] = M.getThermoProp;
    toc
    clear M;
    %filename = sprintf('CoherentDataNEW3_G_Global_Step_wp%.0f_ws%.0f_g%.2f_kh%.4f_kc%.4f.mat',wp,ws,g,kh,kc);
    %save(filename);
end

%% incoherent run

dimpI = 55;

for i = 1:NI
    tic
    fprintf('Incoherent, step %i -- ',i);
    M = modelIncoherentEngine(dimpI,ws,wp,g,km(i),kh,kc,Th,Tc,dissType);
    M.findSS;
    [QmI(i), QhI(i), QcI(i), QbaI(i)] = M.getThermoProp;
    WI(i) = QmI(i) - QbaI(i); %negative when output
    toc
    clear M;
end


%% save the shit
filename = sprintf('cohVSincoh_rates_wp%.0f_ws%.0f_g%.2f_kh%.4f_kc%.4f_nc%.2f.mat',wp,ws,g,kh,kc,nc);
save(filename);


