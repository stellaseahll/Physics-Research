dB = 50;
J = (dB-1)/2;
m = (-J:J).';
rate = sqrt((J-m).*(J+m+1));
rate = rate(1:end-1);
gt = pi/2./rate;
s = chargeBattery(-1,2,dB,1,3);
s.runDiffInteraction(gt);
Es = s.getEnergy();
E2s = s.getE2();
p = chargeBattery(-1,2,dB,1,2);
p.runSameInteraction(10*length(gt),pi/2);
Ep = p.getEnergy();
E2p = p.getE2();