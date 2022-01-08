ppenclear; clc;
delta = logspace(-3,-0.5,51);
g = 1e-3;
E1 = 1;
E2 = 1+delta;
T1 = 1;
T2 = 2;
k1 = 0.01;
k2 = 0.01;
Btype = 'o';
for i = 1:length(g)
    for j = 1:length(E2)
%         fprintf('(%d,%d)\n',i,j);
        Dl=twoQubitDiode(g(i),[T1 T2],[k1 k2],[1 E2(j)],Btype,'l');
        Dg=twoQubitDiode(g(i),[T1 T2],[k1 k2],[1 E2(j)],Btype,'g');
        Dp=twoQubitDiode(g(i),[T1 T2],[k1 k2],[1 E2(j)],Btype,'p');
        Ql(i,j) = Dl.Q(2);
        Qg(i,j) = Dg.Q(2);
        Qp(i,j) = Dp.Q(2);
       
    end
end

% title(sprintf('E1=%.2f,E2=%.2f,T1=%d,T2=%d',E1,E2,T1,T2));
        loglog(delta,Ql,'o');hold on;
       loglog(delta,Qp,'^');
    loglog(delta,Qg,'*'); 
    legend('local','coarse grained','global');
    xlabel('detuning'); 
ylabel('heat flow');
%     
%     xlabel('coupling strength g');
% ylabel('heat flow');
% title(sprintf('E1=%.2f,E2=%.2f,T1=%d,T2=%d',E1,E2,T1,T2));
%         loglog(g,Ql,'o'); hold on;
%        loglog(g,Qp,'^');
%     loglog(g,Qg,'*'); hold off;
%     legend('local','coarse grained','global');