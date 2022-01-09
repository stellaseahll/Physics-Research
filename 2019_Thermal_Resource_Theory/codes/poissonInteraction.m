clc;clf;clear;
%% Initialisation
sx = [0 1; 1 0];
sp = [0 0; 1 0];
sm = [0 1; 0 0];
sz = [-1 0; 0 1];
g = 0.1;
ws = 1;
wb = 1;
T = 1;
Nsteps = 1000;
gamma = 1;
% Distribution: gamma*exp(-gamma*tWait), tWait: waiting time
% Cumulative distribution: 1- exp(-gamma*tWait)
p = rand(1,Nsteps); 
tWait = log(1-p)/(-gamma);
exptp = exp(1i*ws*tWait);
exptm = exp(-1i*ws*tWait);
rhoS = [0 0; 0 1];
rhoB = [1 0; 0 exp(-wb/T)];
rhoB = rhoB/trace(rhoB);
%% Interaction Hamiltonian
% Assume qubit 1 and 2 on resonant, only need to consider Hint in
% interaction picture:

Hint = g*(kron(sp,sm)+kron(sm,sp));
[a b] = eig(Hint);
Uint = zeros(4);
for j = 1:4
    Uint = Uint + a(:,j)*a(:,j)'*exp(-1i*b(j,j));
end


for j = 1:Nsteps
    p0(j) = rhoS(1,1);
    p01(j) = rhoS(1,2);
    rhoS = ptrace(Uint*kron(rhoS,rhoB)*Uint',2,2);
end
plot(cumsum([0 tWait(1:end-1)]),p0)
