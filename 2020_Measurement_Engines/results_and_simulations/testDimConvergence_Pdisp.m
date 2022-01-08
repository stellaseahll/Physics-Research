% get steady state for engine
%test new function where P operator is approximated by displacements
%seems to make no difference, it's just a factor 2 faster

dimp = 20:60; %dimension of oscillator
nh = 1.0; %hot bath excitation
nc = 0.01; %cold excitation
kh = 1e-4; %hot bath rate
kc = 0.001; %cold bath rate
km = 2e-4; %measurement rate
wp = 0.02; %oscillator frequency
ws = 1.0; %spin frequency
g = 0.01; %spin-oscillator coupling

rhoSS = cell(length(dimp),1);
W = zeros(length(dimp),1);
Qh = zeros(length(dimp),1);
Qc = zeros(length(dimp),1);
for j = 1:length(dimp)
    tic
    fprintf(' dim = %i -- ',dimp(j))    
    [Wtmp,Qhtmp,Qctmp,sprhoSS] = getSteadyEngineXbasis_Pdisp(dimp(j),nh,nc,kh,kc,km,wp,ws,g);
    W(j) = Wtmp;
    Qh(j) = Qhtmp;
    Qc(j) = Qctmp;
    rhoSS{j} = full(sprhoSS);
    toc
end

save('testDim_Pdisp.mat');