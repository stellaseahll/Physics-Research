% get steady state for engine

dimp = 81; %dimension of oscillator
nh = 1.0; %hot bath excitation
nc = 0.01; %cold excitation
kh = logspace(-5,-3,50); %hot bath rate
kc = 0.001; %cold bath rate
km = logspace(-5,-3,50); %measurement rate
wp = 0.02; %oscillator frequency
ws = 1.0; %spin frequency
g = 0.01; %spin-oscillator coupling

sigma = 0.1 %sharp POVM, should give abt same as projective
z0 = 0;

rhoSS = cell(length(kh),length(km));
W = zeros(length(kh),length(km));
Qh = zeros(length(kh),length(km));
Qc = zeros(length(kh),length(km));
for j = 1:length(kh)
    for k = 1:length(km)
        [Wtmp,Qhtmp,Qctmp,sprhoSS] = getSteadyEngineXbasis_POVM(dimp,nh,nc,kh(j),kc,km(k),wp,ws,g,z0,sigma);
        W(j,k) = Wtmp;
        Qh(j,k) = Qhtmp;
        Qc(j,k) = Qctmp;
        rhoSS{j,k} = full(sprhoSS);
    end
end

save('run1_POVM.mat');