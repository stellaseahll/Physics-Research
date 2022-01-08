% To test thermalization of qubit (sz) + HO(a,adg) without measurement
% Differs from simply using (a,adg) as dissipators to cold bath, instead
% use b = idQ\otimes a + g/w sz\otimes idHO. 
% Intuitively, this should lead to rhoss = pg |gXg| D(g/2w)rhoth D'(g/2w) +
% pe |eXe| D(-g/2w)rhoth D'(-g/2w) when ws >> g. Otherwise, ws should be
% affected by the dissipator b as well. 
% Factor of 1/2 in displacement op because <sz> = [-1/2,1/2]
%
% Have checked the following: 
% 1) If steady state of qubit is |g> or |e>, harmonic oscillator is in the
% correct displaced thermal state
% 2) If g=0, get product thermal states
% 3) g<<ws, the steady state of qubit = thermal state. Problem is with the
% harmonic oscillator. Tried increasing the dimensions but that didnt seem
% the issue so along as alpha is not too big.
% 4) somehow changing kc affects the steady state also, can't see why that's
% the case
% 5) rhoSS does not even commute with H, but rhoExpect does
% 1) and 2) apply to both ways of computing displacement op (see code)


%clear;clc;
%% Dimensions
dp = 20;

%% Rates and Parameters
ws = 1;
wp = 0.01;
g = 0.01;
Th = 1.0;
Tc = 0.001;
nh = 1/(exp(ws/Th)-1);
nc = 1/(exp(wp/Tc)-1);
kh = 0.0041;
kc = 0.01;
alpha = g/wp/2;

%% Operators for HO
N = diag(0:dp-1);
a = diag(sqrt(1:dp-1),1);
adg = a';
idHO = eye(dp);
% D = getDispOp(dp,alpha); %slower
D=expm(alpha*adg-alpha*a); %faster but errors at high dim

%% Operators for qubit
sz = [-1 0; 0 1]/2;
sp = [0 0; 1 0];
sm = sp';
sx = sp+sm;
idQ = eye(2);


%% Operator for the resolved jump op for bath
b = g/wp*kron(sz,idHO)+kron(idQ,a); 
bdg = b';

%% Thermal States for Comparisons
% Qubit
rhoQ = expm(-ws*sz/Th);
rhoQ = rhoQ/trace(rhoQ);
% Harmonic oscillator
rhoTh = expm(-N*wp/Tc);
rhoTh = rhoTh/trace(rhoTh);
rhoThP = D*rhoTh*D'; %Displaced thermal state in +x if qubit in ground 
rhoThM = D'*rhoTh*D; %Displaced thermal state in -x if qubit in excited
%Combined
rhoExpect = (nh+1)/(2*nh+1)*kron([1 0; 0 0],rhoThP) + nh/(2*nh+1)*kron([0 0; 0 1],rhoThM);
%% Hamiltonian and Dissipators
H = wp*kron(idQ,N) + g*kron(sz,a+adg) + ws*kron(sz,idHO);
U = -1i*(leftMultiply(H) -rightMultiply(H));
Lh =  kh*(nh+1)*(lrMultiply(kron(sm,idHO)) - 0.5*leftMultiply(kron(sp*sm,idHO)) - 0.5*rightMultiply(kron(sp*sm,idHO))) ...
            + kh*(nh)*(lrMultiply(kron(sp,idHO)) - 0.5*leftMultiply(kron(sm*sp,idHO)) - 0.5*rightMultiply(kron(sm*sp,idHO)));
Lc = kc*(nc+1)*(lrMultiply(b) - 0.5*leftMultiply(bdg*b) - 0.5*rightMultiply(bdg*b)) ...
            + kc*nc*(lrMultiply(bdg) - 0.5*leftMultiply(b*bdg) - 0.5*rightMultiply(b*bdg));
Ltot = Lh+Lc+U;

%% Get steady state
rhoSS = reshape(null(Ltot),2*dp,2*dp);
rhoSS = rhoSS/trace(rhoSS);
rhoSSg = rhoSS(1:dp,1:dp)/trace(rhoSS(1:dp,1:dp)); %state conditioned on |g>, should expect rhoSSg = rhoThPp
rhoSSe = rhoSS(dp+1:end,dp+1:end)/trace(rhoSS(dp+1:end,dp+1:end)); %state conditioned on |e>, should expect rhoSSe = rhoThPm
[rhoSSQ,rhoSSHO] = ptrace(rhoSS,2,dp);

