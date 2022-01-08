load('F2randpts.mat');
figure;
Fth = -1./(1+2*N).^2./(1+N)./N;
f = F./Fth/2;
idx = f>1;
f = f(idx);
Gtmp = G(idx);
Ntmp = N(idx);
[f sortidx] = sort(f);
Gtmp = Gtmp(sortidx);
Ntmp = Ntmp(sortidx);
subplot(2,2,1); scatter(Gtmp,Ntmp,1,f);
title('colored region has F2/2Fth > 1')
xlabel('gamma'),ylabel('nbar')
colorbar;

f = F./F0;
idx = f>1;
f = f(idx);
Gtmp = G(idx);
Ntmp = N(idx);
[f sortidx] = sort(f);
Gtmp = Gtmp(sortidx);
Ntmp = Ntmp(sortidx);
subplot(2,2,2); scatter(Gtmp,Ntmp,1,f);
title('colored region has F2/F2_{ground} > 1')
xlabel('gamma'),ylabel('nbar')
colorbar;

Fth = -1./(1+2*N).^2./(1+N)./N;
f = F./Fth/2;
idx = f<1;
idx = idx(1:100000);
f = f(idx);
Gtmp = G(idx);
Ntmp = N(idx);
[f sortidx] = sort(f);
Gtmp = Gtmp(sortidx);
Ntmp = Ntmp(sortidx);
subplot(2,2,3); scatter(Gtmp,Ntmp,1,f);
title('colored region has F2/2Fth < 1')
xlabel('gamma'),ylabel('nbar')
colorbar;

f = F./F0;
idx = f<1;
idx = idx(1:100000);
f = f(idx);
Gtmp = G(idx);
Ntmp = N(idx);
[f sortidx] = sort(f);
Gtmp = Gtmp(sortidx);
Ntmp = Ntmp(sortidx);
subplot(2,2,4); scatter(Gtmp,Ntmp,1,f);
title('colored region has F2/F2_{ground} < 1')
xlabel('gamma'),ylabel('nbar')
colorbar;