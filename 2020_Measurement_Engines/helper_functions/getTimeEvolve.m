%measurement feedback with 2 qubits
clear;clc;
%% Defining Operators
sp = [0 0; 1 0];
sm = sp';
sz = [-1 0; 0 1];
sx = sp+sm;
sy = [0 1i; -1i 0];
s0 = sm*sp;
s1 = sp*sm;
id = eye(2);

sp1 = kron(sp, id);
sm1 = kron(sm, id);
sz1 = kron(sz, id);
sx1 = kron(sx, id);
s01 = kron(s0, id);
s11 = kron(s1, id);
sy1 = kron(sy, id);
sp2 = kron(id, sp);
sm2 = kron(id, sm);
sz2 = kron(id, sz);
sx2 = kron(id, sx);
sy2 = kron(id, sy);
s02 = kron(id, s0);
s12 = kron(id, s1);

%% Rates
gamma = 0.01; %measurement+feedback
kh = 0.001; %hot bath rate
kc = 0.01; %cold bath rate
g = 0.01; %interaction coupling
% tau = 0.001;
% k = logspace(-2,0,100);

%% System and Probe Specs
wh = 1;
wc = 0.1;
Th = 1;
Tc = 0.01;
nh = 1/(exp(wh/Th)-1);
nc = 1/(exp(wc/Tc)-1);
Lh =  kh*(nh+1)*(lrMultiply(sm1) - 0.5*leftMultiply(sp1*sm1) - 0.5*rightMultiply(sp1*sm1)) ...
            + kh*(nh)*(lrMultiply(sp1) - 0.5*leftMultiply(sm1*sp1) - 0.5*rightMultiply(sm1*sp1));
Lc1 = (lrMultiply(sm2) - 0.5*leftMultiply(sp2*sm2) - 0.5*rightMultiply(sp2*sm2));
Lc2 = (lrMultiply(sp2) - 0.5*leftMultiply(sm2*sp2) - 0.5*rightMultiply(sm2*sp2));
Lc = kc*(nc+1)* Lc1 + kc*(nc)*Lc2;
% % %Conditional unitaries
M0 = s02; %measurement of ground state of probe and do nothing
M1 = sx1*s12; %measurement of excited state of probe and rotate system
Lm = gamma*(lrMultiply(M0)-0.5*leftMultiply(M0'*M0)-0.5*rightMultiply(M0'*M0) +...
    lrMultiply(M1)-0.5*leftMultiply(M1'*M1)-0.5*rightMultiply(M1'*M1));

Hint = g*(s11*sx2);
Uint = -1i*(leftMultiply(Hint) -rightMultiply(Hint));
H0 = wh/2*sz1 + wc/2*sz2;
Ufree = -1i*(leftMultiply(H0) -rightMultiply(H0));
L = Lh+Lc+Lm+Ufree+Uint;

t = (0:0.01:10)*(2*pi*100);
rho0 = kron([1 0; 0 0]/2,[1 0; 0 0]/2);
rho0 = reshape(rho0,4*4,1);
for j = 1:length(t)
    rhoSt(:,:,j) = reshape(expm(L*t(j))*rho0,4,4);
    expZ1(j) = trace(rhoSt(:,:,j)*sz1);
    expZ2(j) = trace(rhoSt(:,:,j)*sz2);
    expX1(j) = trace(rhoSt(:,:,j)*sx1);
    expX2(j) = trace(rhoSt(:,:,j)*sx2);
    QH(j) = trace((kh*(nh+1)* diss(rhoSt(:,:,j),sm1) + kh*(nh)* diss(rhoSt(:,:,j),sp1))*(H0+Hint));
    QC(j) = trace((kc*(nc+1)* diss(rhoSt(:,:,j),sm2) + kc*(nc)* diss(rhoSt(:,:,j),sp2))*(H0+Hint));
    W(j) = trace((gamma*diss(rhoSt(:,:,j),M0))*(H0+Hint))+trace((gamma*diss(rhoSt(:,:,j),M1))*(H0+Hint));
end