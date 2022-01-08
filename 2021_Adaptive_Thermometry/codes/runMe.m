Tmin = 1; 
Tmax = 10;
Np = [2 4 8 16 20 30 40 45];
N = 1e4; %total number of qubits 

%% GET DATA
%oldBound = Eq.24, newBound = Mathias' bound, n = total qubits (x axis)

for i = 1:length(Np)
       [n,oldBound,newBound,oldE(i),newE(i)] = getBound(Tmin,Tmax,Np(i),N);
       loglog(n,newBound); hold on;
%         subplot(length(Np),2,2*i-1);
%         loglog(n,newBound,'color','k','linewidth',2); hold on;
%         loglog(n,oldBound,'color','r','linewidth',2);
%         legend('new bound','old bound');
%         xlabel('no. of qubits');
%         ylabel('mle');
%         title(sprintf('interacting qubits = %d',Np(i)));
%         axis([1 1e4 1e-5 1e-1])
%         subplot(length(Np),2,2*i);
%         semilogx(n,oldBound-newBound,'color','k','linewidth',2); hold on;
%         plot([1 1e4],[0 0],'k--','linewidth',0.5);
%         axis tight
%         xlabel('no. of qubits');
%         ylabel('difference mle');
%         title(sprintf('interacting qubits = %d',Np(i)));
end

% It seems like the old bound is tighter since all newBound<oldBound for
% Np=4,8.

