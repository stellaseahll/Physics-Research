% get steady state for engine

dimp = 61; %dimension of oscillator
nh = 1.0; %hot bath excitation
nc = 0.01; %cold excitation
kh = logspace(-5,-3,21); %hot bath rate
kc = 0.001; %cold bath rate
km = logspace(-5,-3,21); %measurement rate
wp = 0.02; %oscillator frequency
ws = 1.0; %spin frequency
g = 0.01; %spin-oscillator coupling

rhoSS = cell(length(kh),length(km));
W = zeros(length(kh),length(km));
Qh = zeros(length(kh),length(km));
Qc = zeros(length(kh),length(km));
for j = length(kh):-1:1
    j
    for k = 1:length(km)
        k
        [Wtmp,Qhtmp,Qctmp,sprhoSS] = getSteadyEngineXbasis(dimp,nh,nc,kh(j),kc,km(k),wp,ws,g);
        W(j,k) = Wtmp;
        Qh(j,k) = Qhtmp;
        Qc(j,k) = Qctmp;
        rhoSS{j,k} = full(sprhoSS);
        save('dim61part2.mat');
    end
end

