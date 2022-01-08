%https://arxiv.org/pdf/1711.03299.pdf
clear;clc;clf;
%% Global
g = logspace(-3,0,50);
g(end) = [];
chi = 0.001;
Es = [5 1 4];
Ts = [2 1 8];
ks = chi*Es;
for j = 1:length(g)
m = model3QubitsFridgeOld(g(j),Ts,ks,Es,'o','x','g');
m.findSS;
Eg(j) = m.getEntanglement();
Cg(j,:) = m.getAllCoherence();
QHg(j) = m.getTotalHeatFlow(1);
QCg(j) = m.getTotalHeatFlow(2);
QWg(j) = m.getTotalHeatFlow(3);
end

%% Partial
for j = 1:length(g)
m = model3QubitsFridgeOld(g(j),Ts,ks,Es,'o','x','p');
m.findSS;
Ep(j) = m.getEntanglement();
QHp(j) = m.getTotalHeatFlow(1);
QCp(j) = m.getTotalHeatFlow(2);
QWp(j) = m.getTotalHeatFlow(3);
Cp(j,:) = m.getAllCoherence();
end

%% Local
for j = 1:length(g)
m = model3QubitsFridgeOld(g(j),Ts,ks,Es,'o','x','l');
m.findSS;
El(j) = m.getEntanglement();
Cl(j,:) = m.getAllCoherence();
QHl(j) = m.getTotalHeatFlow(1);
QCl(j) = m.getTotalHeatFlow(2);
QWl(j) = m.getTotalHeatFlow(3);
end

semilogx(g,El,g,Eg,g,Ep);
legend('local','global','partial');