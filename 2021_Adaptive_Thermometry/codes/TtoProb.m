function p = TtoProb(T,d,E)

p = 1./(1+d*exp(-E./T));
