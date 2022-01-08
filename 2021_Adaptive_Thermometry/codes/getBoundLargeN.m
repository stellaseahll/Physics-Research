function [n,newBound,newE] = getBoundLargeN(Tmin,Tmax,Np,N)
%Script assumes 2 level spectrum after interaction since it can be shown
%numerically that the optimal is still 1 ground state + Nd degenerate states
%where Nd = 2^Np-1, Np = number of INTERACTING qubits, N is the total
%number of qubits

d = 2^Np-1;
NT = 400; %number of T intervals for numerical integration (not so impt here)
alpha = 1; %lets fix the prior
n = Np:Np:N;

[F C1 newE] = runParametersLargeN(Tmin,Tmax,NT,alpha,d);
newBound = 1./(F+(n/Np)*C1);
