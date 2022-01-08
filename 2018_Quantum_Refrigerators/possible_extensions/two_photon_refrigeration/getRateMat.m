function p = getRateMat(g,Th,Tc,Tw,kh,kc,kw,Eh,Ec,Ew)

E = [Eh Ec Ew];
T = [Th Tc Tw];
k = [kh kc kw];
Ng = model3SpinsFridgeTwoPhotonWC(g,T,k,E,'o','r','g',1);
Ng.getRate;
p = full(Ng.W);