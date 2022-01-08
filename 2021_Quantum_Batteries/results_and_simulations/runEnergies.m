clear;
q = 0.25;
theta = pi/2;
tend = 1000;
%bounded
s = chargeBatteryNEW(q,2,201,theta/2,2,'incoh');
t = chargeBatteryNEW(q,2,201,theta/2,2,'coh');
s.runSameInteraction(tend);
% sDist = s.Edist();
sE2 = s.getE2();
sE = s.getEnergy();
sErg = s.getErgotropy();
sDeErg = s.getDephasedErgotropy();
sSysErg = s.getSystemErgotropy();
sEff = sErg/sSysErg;
sPure = s.getPurity();
t.runSameInteraction(tend);
% tDist = t.Edist();
tE2 = t.getE2();
tE = t.getEnergy();
tErg = t.getErgotropy();
tDeErg = t.getDephasedErgotropy();
tSysErg = t.getSystemErgotropy();
tEff = tErg/tSysErg;
tPure = t.getPurity();
clear s
clear t
save('fig3_bounded.mat');
%unbounded
s = chargeBatteryNEW(q,2,201,theta/2,1,'incoh');
t = chargeBatteryNEW(q,2,201,theta/2,1,'coh');
s.runSameInteraction(tend);
% sDist = s.Edist();
sE2 = s.getE2();
sE = s.getEnergy();
sErg = s.getErgotropy();
sDeErg = s.getDephasedErgotropy();
sSysErg = s.getSystemErgotropy();
sEff = sErg/sSysErg;
sPure = s.getPurity();
t.runSameInteraction(tend);
% tDist = t.Edist();
tE2 = t.getE2();
tE = t.getEnergy();
tErg = t.getErgotropy();
tDeErg = t.getDephasedErgotropy();
tSysErg = t.getSystemErgotropy();
tEff = tErg/tSysErg;
tPure = t.getPurity();
clear s
clear t
save('fig3_unbounded.mat');
