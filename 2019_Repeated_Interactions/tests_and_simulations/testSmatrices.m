%run script to get data to compare scattering matrix with numerics
%paper plot in paperSimCompare
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
omegaS = 20;
dimB = 2;
delta = 0;
omegaB = omegaS-delta;
TB = 15;
gsb = diag([1,0.5,0.2,0,0]);
gamma = 50/1000;
dtStep = pi/gamma/100;
Nstep = 800;
tEnd = Nstep*dtStep;
Ntrials = 1e3;
tau = [0.05 0.1 1]/gamma;
t = 0:dtStep:tEnd;
rhoX = eye(dimS)/dimS;
rhoTh = exp(-(1:dimS)*omegaB/TB);
rhoTh = diag(rhoTh)/sum(rhoTh);
rho2 = exp(-(1:dimS)*omegaS/TB);
rho2 = diag(rho2)/sum(rho2);
rho0 = 1/2*[1 1; 1 1];
c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,0.1,gsb,rho0);
c.findSSwithoutSim();
clf;
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
C = colormap(lines);
% figure;
for j = 1:length(tau)
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,tau(j),gsb,rho0);
    c.findSSwithoutSim();
    [rhoSt{j}, obst{j}] = c.runTheWorkers(Ntrials,dtStep,Nstep,obs);
    [rhoStS{j}, obstS{j}] = c.runMEScatter(t,obs);
    [rhoStEik{j}, obstEik{j}] = c.runMESeik(t,obs);
    [rhoStUint{j}, obstUint{j}] = c.runMEApprox(t,obs);
%     save('simResult.mat');
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
    subplot(2,1,1);hold on;
    plot(t,obstS{j}(1,:),'color',C(j,:),'linewidth',3); plot(t,obst{j}(1,:),'color',C(j,:),'linestyle','--'); hold on;
    subplot(2,1,2);hold on;
    plot(t,abs(obstS{j}(3,:)),'color',C(j,:),'linewidth',3); plot(t,abs(obst{j}(3,:)),'color',C(j,:),'linestyle','--'); hold on;
end
% figure; subplot(2,2,1); hold on;
% plot(t,obstS{1}(1,:),'k',t,obstEik{1}(1,:),'b',t(1:100:end),obstUint{1}(1,1:100:end),'ro');
% subplot(2,2,2); hold on;
% plot(t,abs(obstS{1}(3,:)),'k',t,abs(obstEik{1}(3,:)),'b',t,abs(obstUint{1}(3,:)),'r');
% subplot(2,2,3); hold on;
% plot(t,obstS{2}(1,:),'k',t,obstEik{2}(1,:),'b',t(1:100:end),obstUint{2}(1,1:100:end),'ro');
% subplot(2,2,4); hold on;
% plot(t,abs(obstS{2}(3,:)),'k',t,abs(obstEik{2}(3,:)),'b',t,abs(obstUint{2}(3,:)),'r');