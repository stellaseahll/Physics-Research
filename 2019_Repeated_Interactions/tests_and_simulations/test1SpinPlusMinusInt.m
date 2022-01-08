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
delta = 1;
omegaB = omegaS-delta;
TB = 15;
gsb = diag([0.5,0.3,0.2,0,0]);
gamma = 1e-3;
dtStep = 0.01/gamma;
tEnd = 10/gamma;
Nstep = tEnd/dtStep;
Ntrials = 10000000;
tau = [1e1 1e2 1e3];
t = 0:dtStep:tEnd;
rhoX = eye(dimS)/dimS;
rhoTh = exp(-(1:dimS)*omegaB/TB);
rhoTh = diag(rhoTh)/sum(rhoTh);
rho2 = exp(-(1:dimS)*omegaS/TB);
rho2 = diag(rho2)/sum(rho2);
rho0 = 1/2*[1 1; 1 1];
c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,0.1,gsb,rho0);
c.findSSwithoutSim();

% H = c.H0{1} + c.Hint{1};
% E = unique(abs(eig(H)));
% rat = E(2)/E(1);
% dtRes = 2*pi/E(1);
% step = (0.1:0.02:4);
% dtInt = step*dtRes;
% %Analytical results
% g = gpm*2;
% lambda = sqrt(delta^2 + g^2);
% A = delta + lambda;
% N = g^2 + A^2;
% k = 2*g*A/N*sin(lambda*dtInt/2);
% C = exp(-1i*(delta+lambda)*dtInt/2)*g^2/N + exp(-1i*(delta-lambda)*dtInt/2)*A^2/N;
obs{1} = [-0.5 0; 0 0.5];
obs{2} = [0 1; 1 0];
obs{3} = [0 1; 0 0];
for j = 1:length(tau)
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,tau(j),gsb,rho0);
    [rhoSt{j}, obst{j}] = c.runTheWorkers(Ntrials,dtStep,Nstep,obs);
    [rhoStME{j}, obstME{j}] = c.runMEScatter(t,obs);
    save('simResult.mat');
%     plot(t,obst{j}(1,:),t,obstME{j}(1,:)); hold on;
%     plot(t,abs(obst{j}(3,:)),t,abs(obstME{j}(3,:))); hold on;
%     clear c
%     c.getSteadyVal();
%     L1 = 0.25*(1+abs(C)^2)*(eye(4)+c.rightMultiply(sz)*c.leftMultiply(sz)) +...
%         0.25*(p1-p0)*(1-abs(C)^2)*(c.rightMultiply(sz)+c.leftMultiply(sz)) +...
%         0.25*(C+C')*(eye(4) -c.rightMultiply(sz)*c.leftMultiply(sz)) +...
%         0.25*(C-C')*(c.rightMultiply(sz)-c.leftMultiply(sz)) + ...
%         k^2*p0*(c.rightMultiply(sp)*c.leftMultiply(sm)) +...
%         k^2*p1*(c.rightMultiply(sm)*c.leftMultiply(sp));
%     L2 = eye(4);
%     L3 = -1i*omegaS/2*(c.leftMultiply(sz)-c.rightMultiply(sz)) ;
%     L = L3 + gamma*(L1-L2);
%     K1 = k^2*p0*(c.rightMultiply(sp)*c.leftMultiply(sm) - 0.5*c.leftMultiply(sp*sm) - 0.5*c.rightMultiply(sp*sm)) +...
%         k^2*p1*(c.rightMultiply(sm)*c.leftMultiply(sp) - 0.5*c.leftMultiply(sm*sp) - 0.5*c.rightMultiply(sm*sp));
%     K2 = abs(1-C)^2/4*(c.rightMultiply(sz)*c.leftMultiply(sz) - eye(4));
%     K3 = 0.5*1i*imag(C)*(c.leftMultiply(sz)-c.rightMultiply(sz));
%     K = L3 + gamma*(K1+K2-K3);
end
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

% [rhoSt, Jzt, NjumpsAv, t] = c.runTheWorkers(0,0,0); %parallelized
%plot average Jz, compare to thermal bath value
% hold on;
% plot([0,max(t)],[1,1]*sum(c.pthB{1}.*diag(c.JzB)),'--',t,real(Jzt))
% plot([0,max(t)],[1,1]*(rhoS(2,2)- rhoS(1,1))*0.5,'--',t,real(Jzt))
% plot([0,max(t)],[1,1]*(WSS)*avgJump,'--',t,real(Wt))
% WSS*avgJump - mean(Wt(ceil(2/3*Nstep):end))
% plot(t,cumsum(real(Wt)))
% hold on;