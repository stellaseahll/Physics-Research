% clear;clc;
%% Parameters
beta = 1;
ws = 1;
delta = 0.01;
dim = 3;
theta = pi/10;
cos2 = cos(theta)^2;
sin2 = 1-cos2;
isc = 1i*cos(theta)*sin(theta);
alpha = 2;

%% Initial State
rhoi = diag(rand(1,dim));
rhoi = rhoi/trace(rhoi);

% %% Thermal State
H = 0:(dim-1);
% tau = diag(exp(-beta*wp*H));
% Ztau = trace(tau);
% tau = tau/Ztau;
% Etau = wp*sum(H.*diag(tau)');
% Stau = getEntropy(tau,alpha);

% Ht = kron(diag(H),eye(dim))*ws + kron(eye(dim),diag(H))*wp;
% rhot = expm(-beta*Ht);
% Zt = trace(rhot);
% rhot = rhot/Zt;

rhoS = diag(exp(-beta*ws*H));
Zrho = trace(diag(exp(-beta*ws*H)));
rhoS = rhoS/Zrho;
% FUidiv = -log(Ztau)/beta;

Nstep = 100;
Ntrials = 100000;
wp = ws;%randn(Ntrials,Nstep)*delta+ws;
for i = 1:Ntrials
%% Run evolution
rhoi = zeros(3);  rhoi(1) = 1;
%     B = -logm(kron(rhoi,tau));
% rhoI = rhoi;
for t = 1:Nstep
    tau = diag(exp(-beta*(wp+randn*0.1)*H));
    Ztau = trace(tau);
    tau = tau/Ztau;
% %     for alpha = 1:10
% %         Falpha(alpha,t) = getRenyiDiv(rhoi,tauInf)
% %     end
%     Ht = kron(diag(H),eye(dim))*ws + kron(eye(dim),diag(H))*wp;
%     rhot = expm(-beta*Ht);
%     Zt = trace(rhot);
%     rhot = rhot/Zt;
%     B = -logm(kron(rhoi,tau));
    
%     Fsu(i,t) = getRenyiDiv(kron(rhoi,tauInf),rhot,alpha);
%     Fs(i,t) = getRenyiDiv(rhoi,rhoS,alpha);
%     Etau(i,t) = wp*sum(H.*diag(tau)');
%     Stau(i,t) = getEntropy(tau,alpha);
%     FSidiv(i,t) = -log(Zrho)/beta + getRenyiDiv(rhoi,rhoS,alpha)/beta;
%     Erhoi(i,t) = sum(diag(rhoi)'.*H)*ws;
%     Srhoi(i,t) = getEntropy(rhoi,alpha);
%     rhoSUi = kron(rhoi,tau);
%     rhotmp = zeros(dim^2);
%     for j = 1:dim
%         for k = 1:dim
%             rhotmp((j-1)*dim+k,(k-1)*dim+j) = rhotmp((j-1)*dim+k,(k-1)*dim+j)-isc*rhoSUi((j-1)*dim+k,(j-1)*dim+k);
%             rhotmp((k-1)*dim+j,(j-1)*dim+k) = rhotmp((k-1)*dim+j,(j-1)*dim+k)+isc*rhoSUi((j-1)*dim+k,(j-1)*dim+k);
%         end
%     end
%     
%     rhoSUf = cos2*rhoSUi + rhotmp + sin2*kron(tau,rhoi);
%     Fsu2(i,t)  = getRenyiDivND(rhoSUf,rhot,alpha);
%     SSU(t) = getEntropy(rhoSUf,alpha);
%     ESU(t) = trace(Ht*rhoSUf);
%     tauf = cos2*tau + sin2*rhoi;
%     Etauf(t) = wp*sum(H.*diag(tauf)');
%     Stauf(t) = getEntropy(tauf,alpha);
    rhof = cos2*rhoi + sin2*tau;
%     for alpha = 1:4
%         expB(alpha,t) = trace(kron(rhof,tauf)*B^(alpha))/norm(B^(alpha))/alpha;
%     end
%     Erhof(t) = ws*sum(H.*diag(rhof)');
%     Srhof(t) = getEntropy(rhof,alpha);
%     Idiv(t) = getRenyiDiv(rhoSUf,kron(rhof,tauf),alpha);
    
    rhoi = rhof;
%     FSfdiv(t) = -log(Zrho)/beta + getRenyiDiv(rhof,rhoS,alpha)/beta;
%     FUfdiv(t) = -log(Ztau)/beta + getRenyiDiv(tauf,tau,alpha)/beta;
	Fsuf(i,t) = getRenyiDivND(kron(rhof,tauf),rhot,alpha);
end
% Fdiv(i,:) = FSfdiv;
% I(i,:) = Srhof + Stauf - SSU;
% W(i,:) = ESU - Erhoi - Etau;
% dFS(i,:) = (Erhof-Erhoi) - (Srhof-Srhoi)/beta;
% dFU(i,:) = (Etauf-Etau) - (Stauf-Stau)/beta;
% FSU(i,:) = ESU - SSU/beta;
% FS(i,:) = Erhof - Srhof/beta;
% FU(i,:) = Etauf - Stauf/beta;

end

% hold on; plot(0:Nstep-1,mean(Fsu),'b',0:Nstep-1,mean(Fs),'r--');
% 