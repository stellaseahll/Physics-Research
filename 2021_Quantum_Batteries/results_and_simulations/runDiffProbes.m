q = 0;
theta = pi/4;
dS = 2;
dB = 100;
s = chargeBatteryNEW(0,dS,dB,theta,'ho','incoh');

rates = diag(s.bops{1},-1);
thetas = pi/2./rates;
s.runDiffInteraction(thetas)

t = chargeBatteryNEW(0.5,dS,dB,theta,'ho','coh');
rates = diag(t.bops{1},-1);
thetat = pi/4./rates;
thetat = thetat(1:80);
t.runDiffInteraction(thetat)
% 
% theta = thetat;
% qs = [];
% for i = 1:500
%     fprintf('step %i \n',i);
%     theta = [theta;thetat(end)];
%     t.gt = thetat(end);
%     [q f] = myMinSearch(@(q) optqgiventheta(q,thetat(end),dS,dB,t.rhoB{end}),0,0.5,1e-5);
%     rho = t.createCoherentProbe(q);
%     t.runDiffProbes({rho})
%     qs(i) = q;
% end
% 
% figure();
plot(cumsum([0; thetat]),t.getEnergy)
hold on;
plot(cumsum([0; thetas]),s.getEnergy)
hold off