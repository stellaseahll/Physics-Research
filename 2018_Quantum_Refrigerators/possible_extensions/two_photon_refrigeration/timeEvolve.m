clear;clc;
TC = 1;
TH = 2;
TW = 8;
g = logspace(-3,0,100);
g(end) = [];
kappa = 0.00001;
Es = [5 1 4];
data = cell(1,length(g));
ks = ones(1,3)*kappa.*Es;
tmax = zeros(1,length(g));

for i = 1:length(g)
    tic
    M = model3QubitsFridgeOld(g(i),[TH TC TW],ks,Es,'ohmic','x','p');
    M.findSS();
    rhoIn = M.generateRhoIn;
    y = M.timeEvolve(rhoIn,0.05,1000);
    data{i} = y;
    [a b] = max(y(4,:));
    tmax(i) = y(1,b);
    fprintf('Temp %i -- ',i);
    toc
end
save timeevolve_xxx_partial_EC1_EH5_EW4.mat