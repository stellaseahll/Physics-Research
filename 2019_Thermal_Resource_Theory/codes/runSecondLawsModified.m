clear;clc;
%% Parameters
beta = 1;
ws = 1;
delta = 0.01;
dim = 3;
theta = pi/10;
cos2 = cos(theta)^2;
sin2 = 1-cos2;
isc = 1i*cos(theta)*sin(theta);
alpha = 1;

%% Initial State
rhoi = diag(rand(1,dim));
rhoi = rhoi/trace(rhoi);

% %% Thermal State
H = 0:(dim-1);

rhoS = diag(exp(-beta*ws*H));
Zrho = trace(diag(exp(-beta*ws*H)));
rhoS = rhoS/Zrho;
rhoS = diag(rhoS)';
% FUidiv = -log(Ztau)/beta;

Nstep = 50;
Ntrials = 1;
wp = ws+zeros(Ntrials,Nstep)*0.4;%1.2+randn(Ntrials,Nstep)*0.05;
betap = ones(Ntrials,Nstep)*beta;
rhoi = cell(Nstep,1);
rhoI = zeros(1,3);
rhoI(1) = 1;

for i = 1:Ntrials
  	rhoi{i,1} = rhoI;
end
alpha = 2;
Ss = zeros(Ntrials,1);
Ss(:,1) = getEntropy(diag(rhoI),alpha);
Eb = zeros(Ntrials,1);
dEp=zeros(Ntrials,Nstep-1);
dSs=zeros(Ntrials,Nstep-1);
dB=zeros(Ntrials,Nstep-1);
for i = 1:Ntrials
    for t = 2:Nstep
%         Balpha = (-log(kron(rhoi{i,t-1},tau))).^alpha;
        tau = exp(-betap(i,t)*wp(i,t)*H);
        tau = tau/sum(tau);
        tauf = cos2*tau + sin2*rhoi{i,t-1};
        rhoi{i,t} = cos2*rhoi{i,t-1} + sin2*tau;
        Ss(i,t) = getEntropy(diag(rhoi{i,t}),alpha);
        dEp(i,t-1) = getEDiv(diag(tauf),betap(i,t),diag(H)*wp(i,t),alpha)-getEDiv(diag(tau),betap(i,t),diag(H)*wp(i,t),alpha);
%         dB(i,t-1) = sum((kron(rhoi{i,t},tauf)-kron(rhoi{i,t-1},tau)).*Balpha);
%         dD(i,t-1) = -getRenyiDiv(rhoi{i,t},tau,1) +  getRenyiDiv(rhoi{i,t-1},tau,1);
%         dEp(i,t-1) = sum((tauf-tau).*H)*wp(i,t);
%         bdEp(i,t-1) = dEp(i,t-1)*betap(i);
%         Ss(i,t) = getEntropy(diag(rhoi{i,t}));
%         bdEs(i,t-1) = betap(i,t)*sum((rhoi{i,t}-rhoi{i,t-1}).*H)*ws;
    end
end

rhoAvg = zeros(Nstep,3);
% 
% for j = 1:Nstep
%   	rhoAvg{j} = zeros(3);
% end

for j = 1:Nstep
    for i = 1:Ntrials
        rhoAvg(j,:) = rhoAvg(j,:) + rhoi{i,j};
    end
    rhoAvg(j,:) = rhoAvg(j,:)/Ntrials;
end
alpha = 1;

for k = 1:length(alpha)
%     Balpha =  -(log(kron(rhoAvg(1,:),rhoS))).^alpha(k)./alpha(k);
%     Balpha(isinf(Balpha)) = 0;
%     Balpha = Balpha/norm(Balpha);
    for  j = 2:Nstep
%         dBavg(k,j-1) = sum(Balpha.*(kron(rhoAvg(j,:),rhoS)-kron(rhoAvg(1,:),rhoS)));
        dSavg(k,j-1) = getEntropy(diag(rhoAvg(j,:)),alpha(k)) - getEntropy(diag(rhoAvg(j-1,:)),alpha(k));
    end
end
% 
% alpha = 1:5:200;
% Dalpha= 0;
% for i = 1:Nstep
%     for j = 1:length(alpha)
%         Dalpha(i,j) = getRenyiDiv(rhoAvg{i},rhoS,alpha(j));
%     end
% end
% % hold on; plot(0:Nstep-1,mean(Fsu),'b',0:Nstep-1,mean(Fs),'r--');
% % 