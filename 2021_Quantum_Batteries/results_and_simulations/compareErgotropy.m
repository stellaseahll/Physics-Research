clear;clc;
negT = -[0.1 1 10];
dS = 2;
for j = 1:length(negT)
    fprintf('%d\n',j)
    s = chargeBattery(negT(j),dS,200,pi/2,1);
    s.runSameInteraction(50);
    E(j,:) = s.getEnergy();
    Erg(j,:) = s.getErgotropy();
    tmp = diff(E(j,:));
    ErgPerQubit(j) = s.getSystemErgotropy();
    dE(j) = tmp(end);

%     tmp = diff(Erg(j,:));
%     if abs((tmp(end)-tmp(end-1))/(tmp(end)+tmp(end-1)))<1e-5 %to check if boundary is reached.
%         dErg(j) = tmp(end);
%     end
end

