clear;
wp = 1;
ws = 100;
% nc = 1;
kc = 0.1;
kh = 0.001;
g = 2.5;
nh = 1;
Th = ws/log(1/nh+1);
Xend =  g/wp/sqrt(1/2/wp);
X = linspace(0.7,Xend,60);
xth = (g/wp)./X;
Tc = wp/2./acoth(2*wp*xth.^2);
nc = 1./(exp(wp./Tc)-1);
dissType = 'g';
dimp = 55;
%for km = [0.01 0.001 0.0001]
for km = [0.1]
    % Tc = logspace(-1,0.5,10);
    for i = 1:length(Tc)
        i
        M = modelIncoherentEngine(dimp,ws,wp,g,km,kh,kc,Th,Tc(i),dissType);
        M.findSS;
        [W(i), Qh(i), Qc(i)] = M.getThermoProp;
        [Wq(i), Qhq(i), Qcq(i)] = M.getThermoPropBare;
    end
    clear M;
    filename = sprintf('IncoherentData_wp%.0f_ws%.0f_g%.2f_km%.4f_kh%.4f_kc%.4f.mat',wp,ws,g,km,kh,kc);
    save(filename);
end