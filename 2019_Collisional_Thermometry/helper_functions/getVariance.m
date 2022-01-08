function [meanB varB] = getVariance(pn)
%given a string of probabilities '0', we compute the variance based on
%measuring the qubits. R=n0/n0+n1 = 1/(1+exp(beta))
%beta = log(1/R -1)
%Problem: beta diverges for n0 = 0... assign R value based on the next
%smallest R/10
%negative temperatures set to infinite temp;

N = 0:(length(pn)-1);
numBit = log2(length(pn));
numOnes = de2bi(N);
numOnes = sum(numOnes');
numZeros = numBit-numOnes;
R = numZeros./numBit;
R(end) = R(end-1)/10;
beta = log(1./R-1)
beta(beta<0) = 0;
meanB = sum(pn.*beta);
meanBSq = sum(pn.*(beta.^2));
varB = meanBSq-meanB^2;