clear;clc;clf;
% Parameters
wS = 1:100;
wB = 1:100;
for n = 1:100
dimS = 2;
omegaS = wS(n);
dimB = 2;
omegaB = wB(n);
TB = 1;
gamma = 1e-4;
gxx = 0.01;
gyy = 0;
gzz = 0; 
gpm = 0;
dtInt = pi/200:pi/200:2*pi;
% c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,0.01,gxx,gyy,gzz,gpm);
% SS = c.findSSwithoutSim();
% [rhoSt, Jzt, NjumpsAv, t] = c.runTheWorkers(10,50,50000);
rhoX = eye(dimS)/dimS;
rhoTh = exp(-(1:dimS)*omegaS/TB);
rhoTh = diag(rhoTh)/sum(rhoTh);
gsb = [gxx gyy gzz gpm];
for j = 1:length(dtInt)
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(j),gsb);
    rhoS = c.findSSwithoutSim();
    c.getSteadyVal();
    trDistX(j)  = norm(rhoS-rhoX);
    trDistTh(j) = norm(rhoS-rhoTh);
    WSS(j) = real(c.WSS);
    
%     USS(j) = real(c.USS);
%     USSs(j) = real(c.USSs);
%     USSb(j) = real(c.USSb);
%     SSS(j) = real(c.SSS);
%     SSSs(j) = real(c.SSSs);
%     SSSb(j) = real(c.SSSb);
end
    idx = find(WSS(2:end-1)<WSS(1:end-2) & WSS(2:end-1)<WSS(3:end))+1;
    T(n) = idx(1)/200;
end

save('oneSpinXXcheckPeriod.mat');
