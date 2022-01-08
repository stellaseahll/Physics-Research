%measurement feedback with 2 qubits
clear;clc;
%% Defining Operators
ds = 2;
dp = 4;
Js = getSpinOp(ds);
Jp = getSpinOp(dp);

sxs = kron(Js{1}, eye(dp));
sys = kron(Js{2}, eye(dp));
szs = kron(Js{3}, eye(dp));
sps = kron(Js{4}, eye(dp));
sms = kron(Js{5}, eye(dp));
Pms = kron(Js{6}, eye(dp));
Pps = kron(Js{7}, eye(dp));

sxp = kron(eye(ds),Jp{1});
syp = kron(eye(ds),Jp{2});
szp = kron(eye(ds),Jp{3});
spp = kron(eye(ds),Jp{4});
smp = kron(eye(ds),Jp{5});
Pmp = kron(eye(ds),Jp{6});
Ppp = kron(eye(ds),Jp{7});

%% Rates
gamma = 0.01; %measurement+feedback
kh = 0.001; %hot bath rate
kc = 0.01; %cold bath rate
g = 0.01; %interaction coupling

%% System and Probe Specs
wh = 1;
wc = logspace(-3,-1,100);
Th = 1;
Tc = logspace(-3,-1,100);
nh = 1/(exp(wh/Th)-1);

Lh =  kh*(nh+1)*(lrMultiply(sms) - 0.5*leftMultiply(sps*sms) - 0.5*rightMultiply(sps*sms)) ...
            + kh*(nh)*(lrMultiply(sps) - 0.5*leftMultiply(sms*sps) - 0.5*rightMultiply(sms*sps));
Lc1 = (lrMultiply(smp) - 0.5*leftMultiply(spp*smp) - 0.5*rightMultiply(spp*smp));
Lc2 = (lrMultiply(spp) - 0.5*leftMultiply(smp*spp) - 0.5*rightMultiply(smp*spp));
% % %Conditional unitaries
M0 = Pmp; %measurement of ground state of probe and do nothing
M1 = sxs*Ppp; %measurement of excited state of probe and rotate system
Lm = gamma*(lrMultiply(M0)-0.5*leftMultiply(M0'*M0)-0.5*rightMultiply(M0'*M0) +...
    lrMultiply(M1)-0.5*leftMultiply(M1'*M1)-0.5*rightMultiply(M1'*M1));
% M = s02 + sx1*s12;
% Lm = gamma*(lrMultiply(M)-0.5*leftMultiply(M'*M)-0.5*rightMultiply(M'*M));

Hint = g*(Pps*sxp);
Uint = -1i*(leftMultiply(Hint) -rightMultiply(Hint));
for i = 1:length(wc)
    for j = 1:length(Tc)
        %% Hamiltonian and Dissipators
        H0 = wh*szs + wc(i)*szp/(dp-1)*2;
        nc = 1/(exp(wc(i)/Tc(j))-1);
        Lc = kc*(nc+1)* Lc1 + kc*(nc)*Lc2;
        Ufree = -1i*(leftMultiply(H0) -rightMultiply(H0));
        L = Lh+Lc+Lm+Ufree+Uint;
        %% Solve Steady
        rho = reshape(null(L),ds*dp,ds*dp);
        rho = rho/trace(rho);
%         [rhoS,rhoP] = ptrace(rho,2,2);
        QH(i,j) = trace((kh*(nh+1)* diss(rho,sms) + kh*(nh)* diss(rho,sps))*(H0+Hint));
        QC(i,j) = trace((kc*(nc+1)* diss(rho,smp) + kc*(nc)* diss(rho,spp))*(H0+Hint));
%         W(i,j) = trace((gamma*diss(rho,M))*(H0+Hint));
        W(i,j) = trace((gamma*diss(rho,M0))*(H0+Hint))+trace((gamma*diss(rho,M1))*(H0+Hint));
    end
end
% figure; contourf(wc,Tc,real(QC)); hold on; colorbar;
% figure;contourf(wc,Tc,real(QH)); hold on; colorbar;
% figure;contourf(wc,Tc,real(W),[min(min(real(W))) 0 max(max(real(W)))]); hold on; colorbar;
figure;contourf(wc,Tc,real(W),[min(min(real(W))) 0 max(max(real(W)))]); hold on; colorbar;

% plot(wc,wc,'linewidth',2);

% clf;plot(real(W)); hold on; plot(real(QH)); plot(real(QC));