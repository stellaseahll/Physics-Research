% 27/9/2018 Maxwell Demon

clear;clc;
gamma = 1; % omega_s/T_H
lb = 500;
la = 200;
chi = linspace(0.001,1,lb); %T_c/T_h
for i = 1:length(chi)
    alpha(:,i) = linspace(chi(i),1,la);%logspace(log10(chi(i)),0,la)';
end %omega_d/omega_s

chi = ones(la,1) * chi;
ps = 1/ (1 + exp(-gamma)) * ones(la,lb); % ground state of system s
pd = 1./ (1 + exp(-gamma*alpha.*(1./chi))); % ground state of demon d (row = fixed alpha)

Wext = (1-ps).*(2*pd-1); % Extracted work from feedback
Wcost = (2*ps-1).*(1-pd) + alpha.*(pd-ps); % Cost for coupling s and d
Wnet = Wext - Wcost;
Qh = pd-ps;
Qc =  alpha.*(pd-ps);
Eff = Wnet./Qh;
% contourf(chi,alpha,Wnet); figure;
% contourf(chi,alpha,Qh)


% semilogx(alpha,Wcost,'k:',alpha,Wext,'k--',alpha,Wnet,'k')
% clf;
% semilogx(chi(Wnet>0),Wnet(Wnet>0)./Qh(Wnet>0),'k'); hold on;
% semilogx(alpha(Wnet>0),Wext(Wnet>0)./Qh(Wnet>0),'k:'); hold on;

% semilogx(alpha(Wnet>0),1-chi*ones(length(alpha(Wnet>0))),'k--');

% for i = 1:length(chi)
%     idx = find(Wnet(:,i)>0);
%     [tmp,idx2] = max(Wnet(idx,i)./Qh(idx,i));
%     eff(i) = tmp;
%     maxP(i) = Wnet(idx2,i);
%     a(i) = alpha(idx(1));
% end
[~, idx] = max(Wnet);

for i = 1:length(chi)
    measEff(i) = pd(idx(i),i); % efficiency of demon = pd
    EffMaxP(i) = Eff(idx(i),i); % efficiency at max power
end
