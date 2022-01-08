clf;
invT = -1;
theta = 0.95*pi;
tend = 1000;
ptheta = sin(theta/2)^2; 
stheta = sin(theta)/2;
v = chargeBattery(invT,2,50,theta/2,1,'coh');
C = abs(v.rhoSth(1,2));
q = v.rhoSth(1,1);
v.runSameInteraction(tend);
t = 0:tend;
% plot(t,s.getEnergy(),'b',t,ptheta*(1-2*q)*t,'ro');
% peak1 = (-ptheta*(2*q-1)+ 2*C*stheta)*t;
% peak2 = (-ptheta*(2*q-1)- 2*C*stheta)*t;
% for i = 1:2:length(t)
%     plot(diag(v.bops{3}),v.Edist(i,:)); hold on;
%     plot(peak1(i)*[1 1],[0,0.1],'k--');
%     plot(peak2(i)*[1 1],[0,0.1],'k--');     axis([v.bops{3}(1,1) v.bops{3}(end,end) 0 0.1]); hold off;
%     Q(i)=getframe();
% 
% end


for i = 1:length(t)
    plot(diag(v.bops{3}),v.Edist(i,:),diag(s.bops{3}),s.Edist(i,:)); hold on;
%     plot(peak1(i)*[1 1],[0,0.1],'k--');
%     plot(peak2(i)*[1 1],[0,0.1],'k--');     axis([v.bops{3}(1,1) v.bops{3}(end,end) 0 0.1]); hold off;
    Q(i)=getframe();
end