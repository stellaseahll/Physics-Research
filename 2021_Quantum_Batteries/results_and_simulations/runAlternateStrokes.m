clear;
q = 0.25;
c = 1;
theta = 0.01*pi;
N = 200;
tmax = 10000;
rho0 = zeros(N,N);
rho0(1) = 1;

s = chargeBattery_mixed(q,2,N,theta/2,4,'coh',rho0);
probeSequence = [];
for j = 1:5
    probeSequence = [probeSequence,ones(1,100),ones(1,100)*2];
end

s.runDiffProbes_cohandincoh(probeSequence);
E1 = s.getEnergy;
Erg1 = s.getErgotropy;
    
t = chargeBattery_mixed(q,2,N,theta/2,4,'coh',rho0);
probeSequence = ones(1,1000);
t.runDiffProbes_cohandincoh(probeSequence);
E2 = t.getEnergy;
Erg2 = t.getErgotropy;
    

