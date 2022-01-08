% run script to get data for different interaction types
% clear;clc;
% Parameters
dimS = 2;
omegaS = 10*(0.8:0.01:1.2);
dimB = 2;
omegaB = 10;
TB = 10;
gamma = 1e-2;
gxx = 1;
gyy = 0;
gzz = 0; 
gpm = 0;
gsb = [gxx, gyy,gzz,gpm];
Ntrials = 100;
dtInt = pi*[0.1:0.01:0.2];
rhoAvg = zeros(dimS,dimS,length(dtInt));
tEnd = 2e5;
numStep = 2000;
dtStep = tEnd/numStep;
avgJump = gamma*dtStep;
% for j = 1:length(dtInt)
%     c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(j),gsb);
%     rhoS(:,:,j) = c.findSSwithoutSim();
%     c.getSteadyVal();
%     WSS(j) = real(c.WSS);
%     USS(j) = real(c.USS);
%     USSs(j) = real(c.USSs);
%     USSb(j) = real(c.USSb);
%     SSS(j) = real(c.SSS);
%     SSSs(j) = real(c.SSSs);
%     SSSb(j) = real(c.SSSb);
%     [rhoSt{j}, Jzt(j,:), Wt(j,:), NjumpsAv(j,:), t] = c.runTheWorkers(100,dtStep,numStep); 
%     for k = 1001:2000
%         rhoAvg(:,:,j) = rhoAvg(:,:,j) + rhoSt{j}(:,:,k)/1000;
%     end
%     WtAvg(j) = mean(Wt(j,1001:2000));
%     traceDist(j) = norm(rhoS(:,:,j)-rhoAvg(:,:,j));
% end

for j = 1:length(omegaS)
    c = collisionModelSpinBath(dimS,omegaS(j),dimB,omegaB,TB,gamma,0.1*pi,gsb);
    rhoS(:,:,j) = c.findSSwithoutSim();
    c.getSteadyVal();
    WSS(j) = real(c.WSS);
    USS(j) = real(c.USS);
    USSs(j) = real(c.USSs);
    USSb(j) = real(c.USSb);
    SSS(j) = real(c.SSS);
    SSSs(j) = real(c.SSSs);
    SSSb(j) = real(c.SSSb);
%     [rhoSt{j}, Jzt(j,:), Wt(j,:), NjumpsAv(j,:), t] = c.runTheWorkers(100,dtStep,numStep); 
%     for k = 1001:2000
%         rhoAvg(:,:,j) = rhoAvg(:,:,j) + rhoSt{j}(:,:,k)/1000;
%     end
%     WtAvg(j) = mean(Wt(j,1001:2000));
%     traceDist(j) = norm(rhoS(:,:,j)-rhoAvg(:,:,j));
end