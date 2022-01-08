clear;
wp = 1; %should be our unit
ws = 100;
% nc = 1;
kc = logspace(-3,-1,100);
km = logspace(-3,-1,101);
kh = 0.001;
nh = 1;
Th = ws/log(1/nh+1);
dissType = 'g';
dimp = 40;

% N = 100;

%% different g, same dim

g = 2.5;
x0 = g/wp/sqrt(2); %value according to paper definition. 
%In the function we use D(g/wp) for D^2, which is, according to the paper,
%D(sqrt(2)*x0)

%2*x0 must be greater than 1.
%don't include xth=1, which amounts to Tc=0
% xth = 1 + (1:N)/N*(2*x0-1); %from 1 (Tc -> 0) to 2*x0
% Tc = wp/2./acoth(xth.^2);
nc = 0.1;
Tc = wp/log(1/nc+1);

for i = 1:length(kc)
    
    for j = 1:length(km)
        fprintf('(%d,%d)\n',i,j);
%         fprintf('km=%1.1e, step %i: xth=%1.2f, nc = %1.2f -- ',km,i,xth(i),nc(i));
        M = modelIncoherentEngine(dimp,ws,wp,g,km(j),kh,kc(i),Th,Tc,dissType);
        M.findSS;
        [Qm(i,j), Qh(i,j), Qc(i,j), Qba(i,j)] = M.getThermoProp;
        W(i,j) = Qm(i,j) - Qba(i,j); %negative when output
        [Wq(i,j), Qhq(i,j), Qcq(i,j)] = M.getThermoPropBare; %just for debug, OBSOLETE
    end

end
clear M;
filename = sprintf('thesisData_wp%.0f_ws%.0f_g%.2f_kh%.4f_nh%.1f_nc%.1f.mat',wp,ws,g,kh,nh,nc);
save(filename);
