%test script single spin coupled to bath with analytical expression
%Load the object with inputs
%Arguments:
% (dimS,omegaS,dimB,omegaB,TB,gamma,dtint,gxx,gyy,gzz,gpm)
%  -> gamma = overall jump rate
%  -> dtint = length of one collision
%  -> gxx,gyy,gzz,gpm coupling rates in interaction Hamiltonian:
%       Hint = gxx*JxS*JxB + gyy*JyS*JyB + gzz*JzS*JzB + gpm*(J-S*J+B + h.c.)
% clear;clc;
% Parameters
dimS = 2;
omegaS = 10;
dimB = 2;
delta = 0;
omegaB = omegaS-delta;
Omega = omegaS + omegaB;
TB = 10;
gamma = 0.01;
gxx = 1;
gyy = 0;
gzz = 0; 
gpm = 0;
gsb = [gxx gyy gzz gpm];
Ntrials = 100;
rho0 = [1 0; 0 0];
obs{1} = [-0.5 0; 0 0.5];
rhoX = eye(dimS)/dimS;
rhoTh = exp(-(1:dimS)*omegaB/TB);
rhoTh = diag(rhoTh)/sum(rhoTh);
rho2 = exp(-(1:dimS)*omegaS/TB);
rho2 = diag(rho2)/sum(rho2);
c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,0.1,gsb,rho0);
c.findSSwithoutSim();

H = c.H0{1} + c.Hint{1};
E = unique(abs(eig(H)));
rat = E(2)/E(1);
dtRes = 2*pi/E(1);
step = (0.1:0.1:4);
% dtInt = logspace(-5,1,100);
% % step = (0.1:0.02:4)
dtInt = pi/E(1)/2;%step*dtRes;
% dtInt = 0.1;
% %Analytical results
g = gxx/2;
lambda = sqrt(delta^2 + g^2);
A = delta + lambda;
N = g^2 + A^2;
k = 2*g*A/N*sin(lambda*dtInt/2);
C = exp(-1i*(delta+lambda)*dtInt/2)*g^2/N + exp(-1i*(delta-lambda)*dtInt/2)*A^2/N;
mu = sqrt(Omega^2 + g^2);
B = Omega + mu;
M = g^2 + B^2;
m = 2*g*B/M*sin(mu*dtInt/2);
D = exp(-1i*(Omega+mu)*dtInt/2)*g^2/M + exp(-1i*(Omega-mu)*dtInt/2)*B^2/M;
sp = [0 0; 1 0];
sm = [0 1; 0 0];
sz = [-1 0; 0 1];
sx = sp+sm;
p0 = c.pthB{1}(1);
p1 = c.pthB{1}(2);
smCoeff = k.^2.*p0 + m.^2.*p1 -k.*m;
spCoeff = k.^2*p1 + m.^2*p0 - k.*m;
xCoeff = k.*m;

for j = 1:length(dtInt)
    Ntrials = 10000;
    tstep = 5;
    numStep = 200;
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(j),gsb,rho0);
    c.prepareSim();
    % [rhoSt, obst, Wt, NjumpsAv,t] = c.runThemBitches(Ntrials,tstep,numStep,obs);
%     [rhoStME, obstME] = c.runMEScatter(t,obs);
%     [rhoStME2, obstME2] = c.runMEApprox(t,obs);
     c.findSS();
%      rhoSscatter = c.rhoSSscatter{1};
     c.getSteadyVal();
     USSs(j) = c.USSs(1);
end
% rhoSscatter = c.rhoSSscatter{1};
% plot(t,real(obst(1,:)),t,real(obstME(1,:)),'ro',t,real(obstME2(1,:)),'k');
% 
% hold on;
% semilogx(real(t),real(obstME(1,:)),real(t),real(obstME2(1,:)),'r'); 
% for j = 1:length(dtInt)
%     c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(j),gsb);
%     c.findSSwithoutSim();
%     c.getLpm;    
%     rhoS(:,:,j) = c.rhoSS{1}; 
%     z1(j) = rhoS(1,1,j);
%     rhoSscatter(:,:,j) = c.rhoSSscatter{1};
%     z2(j) = rhoSscatter(2,2,j);
% %     c.getSteadyVal();
%     trDistX(j)  = norm(rhoS(:,:,j)-rhoX);
%     trDistTh(j) = norm(rhoS(:,:,j)-rhoTh);
%     trDist(j) = norm(rhoS(:,:,j)-rhoSscatter(:,:,j));
%     L1 = 0.25*(abs(D)^2+abs(C)^2)*(eye(4)+c.rightMultiply(sz)*c.leftMultiply(sz)) +...
%         0.25*(p1-p0)*(abs(D)^2-abs(C)^2)*(c.rightMultiply(sz)+c.leftMultiply(sz)) +...
%         0.25*(C*D+(C*D)')*(eye(4) -c.rightMultiply(sz)*c.leftMultiply(sz)) +...
%         0.25*(C*D-(C*D)')*(c.rightMultiply(sz)-c.leftMultiply(sz)) + ...
%         (k^2*p0+m^2*p1)*(c.rightMultiply(sp)*c.leftMultiply(sm)) +...
%         (k^2*p1+m^2*p0)*(c.rightMultiply(sm)*c.leftMultiply(sp)) +...
%         k*m*(c.rightMultiply(sm)*c.leftMultiply(sm)+c.rightMultiply(sp)*c.leftMultiply(sp));
%     L2 = eye(4);
%     L3 = -1i*omegaS/2*(c.leftMultiply(sz)-c.rightMultiply(sz)) ;
%     L = L3 + gamma*(L1-L2);
%     K1 = (k^2*p0+m^2*p1-k*m)*(c.rightMultiply(sp)*c.leftMultiply(sm) - 0.5*c.leftMultiply(sp*sm) - 0.5*c.rightMultiply(sp*sm)) +...
%         (k^2*p1+m^2*p0-k*m)*(c.rightMultiply(sm)*c.leftMultiply(sp) - 0.5*c.leftMultiply(sm*sp) - 0.5*c.rightMultiply(sm*sp));
%     K2 = abs(D-C)^2/4*(c.rightMultiply(sz)*c.leftMultiply(sz) - eye(4)) + k*m*(c.rightMultiply(sx)*c.leftMultiply(sx) - eye(4));
%     K3 = 0.5*1i*imag(C*D)*(c.leftMultiply(sz)-c.rightMultiply(sz));
%     K = L3 + gamma*(K1+K2-K3);
% end
% hold on;
% plot(step,z1,'k',[step(1) step(end)],rhoTh(1)*[1 1],'b--',[step(1) step(end)],rhoTh(2,2)*[1 1],'r--');
% save('multipledtIntLin.mat');
% clf; 
% subplot(3,1,1);
% title('Sb/Ub');hold on;
% plot(dtInt,SSSb./USSb);
% subplot(3,1,2);
% title('Work');hold on;
% plot(dtInt,WSS);
% subplot(3,1,3);
% title('Power');hold on;
% plot(dtInt,WSS./dtInt);
% print('steadyLin.pdf','-dpdf');
% clf; 
% subplot(3,1,1);
% title('Sb/Ub');hold on;
% semilogx(dtInt,SSSb./USSb);
% subplot(3,1,2);
% title('Work');hold on;
% semilogx(dtInt,WSS);
% subplot(3,1,3);
% title('Power'); hold on;
% semilogx(dtInt,WSS./dtInt);
% print('steadyLog.pdf','-dpdf');
% semilogx(dtInt,trDistX,'k','linewidth',3); hold on;
% semilogx(dtInt,trDistTh,'k--','linewidth',3)
% axis([1e-2 1e2 0 0.25]);
% set(gca,'FontSize', 18);
% set(gca,'XTick',[1e-2 1e-1 1e0 1e1 1e2]);
% set(gca,'YTick',0:0.1:0.2);
% xlabel('Interaction Time Interval \deltat_{int}/\omega_S');
% ylabel('|\rho-\sigma|');
% legend('\sigma = \rho^{mix}','\sigma = \rho^{th}')
%run over many trials, either single core or parallelized
% Arguments: (Ntrials, dt_step, Nsteps)
%[rhoSt, Jzt, NjumpsAv, t] = c.runTheBitch(1000,10,1000); %single core

% [rhoSt, Jzt, NjumpsAv, t] = c.runThemBitches(0,0,0); %parallelized
%plot average Jz, compare to thermal bath value
% hold on;
% plot([0,max(t)],[1,1]*sum(c.pthB{1}.*diag(c.JzB)),'--',t,real(Jzt))
% plot([0,max(t)],[1,1]*(rhoS(2,2)- rhoS(1,1))*0.5,'--',t,real(Jzt))
% plot([0,max(t)],[1,1]*(WSS)*avgJump,'--',t,real(Wt))
% WSS*avgJump - mean(Wt(ceil(2/3*Nstep):end))
% plot(t,cumsum(real(Wt)))
% hold on;