clear;
gamma = 1; % omega_s/T_H
chi = 0.05; %T_c/T_h
alpha = 0.5;
ps = 1/ (1 + exp(-gamma)); % ground state of system s
Pd = 1/ (1 + exp(-gamma*alpha*(1/chi))); % ground state of demon d (row = fixed alpha)
Ntrials = 100000000;
i = 1;
pd = Pd;
Wm = 0;
QC = 0;
QH = 0;
W = 0;
while i<=Ntrials
    Wm = Wm + (2*ps-1).*(1-pd) + alpha.*(pd-ps);
    if (rand()<ps) %measure 0
        W = W+0; 
        QH = QH+(pd-ps);
%         QC = 0; %don't put in cold bath
%         pd = 1; % don't use thermal probes
    else
        W = W+(2*pd-1); 
        QH = QH+(pd-ps);
        QC = QC-alpha*pd;
%         pd = Pd; %next use thermal probes
    end
    i = i+1;
end
W = W/Ntrials;
Wm = Wm/Ntrials;
QH = QH/Ntrials;
QC = QC/Ntrials;
Wnet = W-Wm;
Wth = (1-ps)*(2*Pd-1);
QHth = (Pd-ps);
QCth = alpha*(ps-Pd);
Wmth = (2*ps-1).*(1-Pd) + alpha.*(Pd-ps);
Wnetth = Wth-Wmth;

%%%%%%%%%%%%%%%% WITHOUT RECYCLING %%%%%%%%%%%%%%%%%%%%%
% while i<=Ntrials
%     P = rand();
%     if (isThermal)
%         Wm(i) = (2*ps-1).*(1-pd) + alpha.*(pd-ps);
%     else
%         Wm(i) = (1-ps)*alpha;
%     end
%     if (P<ps) %measure 0
%         W(i) = 0; 
%         QH(i) = (pd-ps);
%         QC(i) = alpha*(1-pd);
%     else
%         W(i) = (2*pd-1); 
%         QH(i) = (pd-ps);
%         QC(i) = -alpha*pd;
%     end
%     i = i+1;
% end
