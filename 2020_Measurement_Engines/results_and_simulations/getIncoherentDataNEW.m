clear;
wp = 1; %should be our unit
ws = 100;
% nc = 1;
kc = 0.1;
kh = 0.001;
nh = 1;
Th = ws/log(1/nh+1);
dissType = 'g';
dimp = 55;

N = 100;

%% different g, same dim

g = 2.5;
x0 = g/wp/sqrt(2); %value according to paper definition. 
%In the function we use D(g/wp) for D^2, which is, according to the paper,
%D(sqrt(2)*x0)

%2*x0 must be greater than 1.
%don't include xth=1, which amounts to Tc=0
xth = 1 + (1:N)/N*(2*x0-1); %from 1 (Tc -> 0) to 2*x0
Tc = wp/2./acoth(xth.^2);
nc = 1./(exp(wp./Tc)-1);

for km = [0.1 0.01 0.001 0.0001]
    % Tc = logspace(-1,0.5,10);
    for i = 1:length(Tc)
        tic
        fprintf('km=%1.1e, step %i: xth=%1.2f, nc = %1.2f -- ',km,i,xth(i),nc(i));
        M = modelIncoherentEngine(dimp,ws,wp,g,km,kh,kc,Th,Tc(i),dissType);
        M.findSS;
        [Qm(i), Qh(i), Qc(i), Qba(i)] = M.getThermoProp;
        W(i) = Qm(i) - Qba(i); %negative when output
        [Wq(i), Qhq(i), Qcq(i)] = M.getThermoPropBare; %just for debug, OBSOLETE
        toc
    end
    clear M;
    filename = sprintf('IncoherentDataNEW2_wp%.0f_ws%.0f_g%.2f_km%.4f_kh%.4f_kc%.4f.mat',wp,ws,g,km,kh,kc);
    save(filename);
end

%% different g, same dim

x0 = 2.0;
g = sqrt(2)*wp*x0;

%2*x0 must be greater than 1.
%don't include xth=1, which amounts to Tc=0
xth = 1 + (1:N)/N*(2*x0-1); %from 1 (Tc -> 0) to 2*x0
Tc = wp/2./acoth(xth.^2);
nc = 1./(exp(wp./Tc)-1);

for km = [0.1 0.01 0.001 0.0001]
    % Tc = logspace(-1,0.5,10);
    for i = 1:length(Tc)
        tic
        fprintf('km=%1.1e, step %i: xth=%1.2f, nc = %1.2f -- ',km,i,xth(i),nc(i));
        M = modelIncoherentEngine(dimp,ws,wp,g,km,kh,kc,Th,Tc(i),dissType);
        M.findSS;
        [Qm(i), Qh(i), Qc(i), Qba(i)] = M.getThermoProp;
        W(i) = Qm(i) - Qba(i); %negative when output
        [Wq(i), Qhq(i), Qcq(i)] = M.getThermoPropBare; %just for debug, OBSOLETE
        toc
    end
    clear M;
    filename = sprintf('IncoherentDataNEW2_wp%.0f_ws%.0f_g%.2f_km%.4f_kh%.4f_kc%.4f.mat',wp,ws,g,km,kh,kc);
    save(filename);
end

%% different g, same dim

x0 = 2.5;
g = sqrt(2)*wp*x0;

%2*x0 must be greater than 1.
%don't include xth=1, which amounts to Tc=0
xth = 1 + (1:N)/N*(2*x0-1); %from 1 (Tc -> 0) to 2*x0
Tc = wp/2./acoth(xth.^2);
nc = 1./(exp(wp./Tc)-1);

for km = [0.1 0.01 0.001 0.0001]
    % Tc = logspace(-1,0.5,10);
    for i = 1:length(Tc)
        tic
        fprintf('km=%1.1e, step %i: xth=%1.2f, nc = %1.2f -- ',km,i,xth(i),nc(i));
        M = modelIncoherentEngine(dimp,ws,wp,g,km,kh,kc,Th,Tc(i),dissType);
        M.findSS;
        [Qm(i), Qh(i), Qc(i), Qba(i)] = M.getThermoProp;
        W(i) = Qm(i) - Qba(i); %negative when output
        [Wq(i), Qhq(i), Qcq(i)] = M.getThermoPropBare; %just for debug, OBSOLETE
        toc
    end
    clear M;
    filename = sprintf('IncoherentDataNEW2_wp%.0f_ws%.0f_g%.2f_km%.4f_kh%.4f_kc%.4f.mat',wp,ws,g,km,kh,kc);
    save(filename);
end