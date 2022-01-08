function [sErg, ergTh,q,tend] = runErgTheoretical(theta)

q = 0.1:0.1:0.5;
% theta = pi/2;
ptheta = sin(theta/2)^2;
tend = 200;
t = 0:tend;
for i = 1:length(q)
    ergTh(i,:) = (1-2*q(i))*ptheta*t - 2*sqrt(2*ptheta*t/pi*(1-(1-2*q(i))^2*ptheta));
end
for i = 1:length(q)
    i
    s = chargeBatteryNEW(q(i),2,200,theta/2,1,'incoh');
    s.runSameInteraction(tend);
    sE2(i,:) = s.getE2();
    sE(i,:) = s.getEnergy();
    sErg(i,:) = s.getErgotropy();
    sErg(i,:) = sErg(i,:)-sErg(i,1);
end
