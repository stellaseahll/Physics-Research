clear;clc;
% Parameters
dimS = 2;
omegaS = 10;
dimB = 2;
omegaB = 10;
TB = 10;
gamma = 1e-2;
gxx = 1;
gyy = 1;
gzz = 1; 
gpm = 0;
gsb = [gxx gyy gzz gpm];
Ntrials = 10000;
dtInt = 1;
sigma = [0.01 0.1 1];
%% 
% [rhoSt, Jzt, Wt, NjumpsAv, t] = c.runTheWorkers(100,1000,50000); 
% SS = c.findSSwithoutSim();
% [rhoSt, Jzt, NjumpsAv, t] = c.runTheWorkers(10,50,50000);
numStep = 1000;
rhoSt = cell(length(sigma),1);
St = cell(length(sigma),1);
Et = cell(length(sigma),1);
Jzt = zeros(length(sigma),numStep+1);
Wt = zeros(length(sigma),numStep+1);
t = zeros(length(sigma),numStep+1);
tstep = [2 10 400];
rho0 = [1 0; 0 0];
for j = 1:length(sigma)
    rhoX = eye(dimS)/dimS;
    rhoTh = exp(-(1:dimS)*omegaS/TB);
    rhoTh = diag(rhoTh)/sum(rhoTh);
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt,gsb,rho0);
    [rhoSt{j}, Jzt(j,:), Wt(j,:), St{j}, Et{j}, NjumpsAv(j), t(j,:)] = c.runTheWorkersNoise(Ntrials,tstep(j),numStep,sigma(j));
    save('noisyBath_wS10_wB10_TB10_gjj1_gamma0.01_dtInt1_Ntrials10000.mat');
end
clear;clc;
% Parameters
dimS = 2;
omegaS = 10;
dimB = 2;
omegaB = 10;
TB = 10;
gamma = 1e-2;
gxx = 1;
gyy = 0;
gzz = 0; 
gpm = 0;
gsb = [gxx gyy gzz gpm];
Ntrials = 10000;
dtInt = 1;
sigma = [0.01 0.1 1];
%% 
% [rhoSt, Jzt, Wt, NjumpsAv, t] = c.runTheWorkers(100,1000,50000); 
% SS = c.findSSwithoutSim();
% [rhoSt, Jzt, NjumpsAv, t] = c.runTheWorkers(10,50,50000);
numStep = 1000;
rhoSt = cell(length(sigma),1);
St = cell(length(sigma),1);
Et = cell(length(sigma),1);
Jzt = zeros(length(sigma),numStep+1);
Wt = zeros(length(sigma),numStep+1);
t = zeros(length(sigma),numStep+1);
tstep = [2 10 400];
rho0 = [1 0; 0 0];
for j = 1:length(sigma)
    rhoX = eye(dimS)/dimS;
    rhoTh = exp(-(1:dimS)*omegaS/TB);
    rhoTh = diag(rhoTh)/sum(rhoTh);
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt,gsb,rho0);
    [rhoSt{j}, Jzt(j,:), Wt(j,:), St{j}, Et{j}, NjumpsAv(j), t(j,:)] = c.runTheWorkersNoise(Ntrials,tstep(j),numStep,sigma(j));
    save('noisyBath_wS10_wB10_TB10_gxx1_gamma0.01_dtInt1_Ntrials10000.mat');
end

% 
% clear;clc;
% % Parameters
% dimS = 2;
% omegaS = 10;
% dimB = 2;
% omegaB = 10;
% TB = 10;
% gamma = 1e-2;
% gxx = 1;
% gyy = 1;
% gzz = 1; 
% gpm = 0;
% gsb = [gxx gyy gzz gpm];
% Ntrials = 100;
% dtInt = 0.1;
% sigma = [0.01 0.1 1];
% %% 
% % [rhoSt, Jzt, Wt, NjumpsAv, t] = c.runTheWorkers(100,1000,50000); 
% % SS = c.findSSwithoutSim();
% % [rhoSt, Jzt, NjumpsAv, t] = c.runTheWorkers(10,50,50000);
% numStep = 100;
% rhoSt = cell(length(sigma),1);
% St = cell(length(sigma),1);
% Et = cell(length(sigma),1);
% Jzt = zeros(length(sigma),numStep+1);
% Wt = zeros(length(sigma),numStep+1);
% t = zeros(length(sigma),numStep+1);
% tstep = [200 1000 10000];
% rho0 = [1 0; 0 0];
% for j = 1:length(sigma)
%     rhoX = eye(dimS)/dimS;
%     rhoTh = exp(-(1:dimS)*omegaS/TB);
%     rhoTh = diag(rhoTh)/sum(rhoTh);
%     c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt,gsb,rho0);
%     [rhoSt{j}, Jzt(j,:), Wt(j,:), St{j}, Et{j}, NjumpsAv(j), t(j,:)] = c.runTheWorkersNoise(Ntrials,tstep(j),numStep,sigma(j));
%     save('noisyBath_wS10_wB10_TB10_gjj1_gamma0.01_dtInt1_Ntrials10000.mat');
% end