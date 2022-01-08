%get data for varying g and Tc.
clear;
wp = 1;
ws = 100;
% nc = 1;
kc = 0.1;
km=0.01;
kh = 0.001;
g = linspace(0.5,4,10);
nh = 1;
Th = ws/log(1/nh+1);
% Xend =  g/wp/sqrt(1/2/wp);
% X = linspace(0.7,Xend,30);
% xth = (g/wp)./X;
Tc = 0.01;
nc = 1./(exp(wp./Tc)-1);
dissType = 'g';
dimp = 50;
%for km = [0.01 0.001 0.0001]
for i=1:length(g)
    % Tc = logspace(-1,0.5,10);
    for j = 1:length(Tc)
        fprintf('(%d,%d)\n',i,j);
        M = modelIncoherentEngine(dimp,ws,wp,g(i),km,kh,kc,Th,Tc(j),dissType);
        M.findSS;
        [W(i,j), Qh(i,j), Qc(i,j)] = M.getThermoProp;
        [Wq(i,j), Qhq(i,j), Qcq(i,j)] = M.getThermoPropBare;
    end
    clear M;
end
filename = sprintf('IncoherentData_wp%.0f_ws%.0f_km%.4f_kh%.4f_kc%.4f.mat',wp,ws,km,kh,kc);
save(filename);