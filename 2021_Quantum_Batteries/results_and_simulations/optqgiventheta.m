function f = optqgiventheta(q,theta,dS,dB,rho)

t = chargeBatteryNEW(q,dS,dB,theta,'ho','coh');
t.rhoB{1} = rho;
t.runSameInteraction(1);
tE =  t.getEnergy;
f = tE(1) - tE(2);