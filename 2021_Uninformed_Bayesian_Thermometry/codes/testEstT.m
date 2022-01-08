%test estTprobeGibbs

NT = 300; %number of T vals for T-distributions
NP = 100; %number of qubit probes per experiment
Tmax = 10; %highest T in distribution interval
d = 1;
%T = (1:NT)/NT * Tmax;
T = logspace(-2,1,NT+1);

%choose prior: flat, Jeffreys, ~ T^x
%whichP = 'flat', 'Jeff', x;
whichP = 'Jeff'; %1/T prior

s = estTprobeGibbs(T,NP,d,whichP);

%generates data for NE experiments per underlying true T value
NE = 1000;
Ntrue = 21;
%Ttrue = linspace(0.1,15,Ntrue)+ (0.5*rand(1,Ntrue)-0.25)*0;
Ttrue = 0.1:0.1:10;
s.generateData(NE,Ttrue);

r1 = s.pointEstimate(); T1 = s.lastEstT;
r2 = s.BayesEstimate(); T2 = s.lastEstT;
r3 = s.JesusEstimate(); T3 = s.lastEstT;