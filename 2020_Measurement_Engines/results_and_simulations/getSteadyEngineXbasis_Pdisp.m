function [W,Qh,Qc,sprhoSS] = getSteadyEngineXbasis_Pdisp(dimp,nh,nc,kh,kc,km,wp,ws,g)
%test measurement engine with x projectors, and HO in xbasis
%SN: For damping Lc use displacement operator of smallest position
%displacement to approximate P operator. This might improve convergence,
%but gives wrong results close to cutoff/resolution limit.

%% Operators for qubit
sz = [-1 0; 0 1];
sp = [0 0; 1 0];
sm = sp';
sx = sp+sm;
idQ = eye(2);

%% Operators for HO
%hbar = 1
%mass m=1
%x = (a+adg)/sqrt(2), so unit should be sqrt(hbar/m/wp)
%p = (a-adg)/1i/sqrt(2)
%x interval from -xmax to xmax
%choice critical! if xmax too high, then wrong results for thermal nbar at
%low dimension dp. If too low, then probably wrong results for high nbar...

N = diag(0:dimp-1);
a = diag(sqrt(1:dimp-1),1);
adg = a';
idHO = eye(dimp);
xmax = 10; 
x = linspace(-xmax,xmax,dimp);
dx = 2*xmax/(dimp-1);
p = linspace(-pi/dx,pi/dx,dimp); %fft ordering

%position operator in position repres.
X = (diag(x.'));
%FT basis trafo, from position to momentum repres.
FTM = exp(-1i*p.'*x)/sqrt(dimp);
% %check unitarity
% fprintf( 'FT unitarity: %1.2e, %1.2e \n', norm(FTM'*FTM-eye(dimp)),norm(FTM*FTM'-eye(dimp)) );

%momentum operator in position repres.
P = FTM'*(diag(p.')*FTM);

%shifted position operator due to interaction
Xeff = kron(idQ,X) + g/wp*kron(sz,idHO);

%% Hamiltonians and Dissipators
Hq = ws/2*kron(sz,idHO);
Hp = wp/2*kron(eye(2),P*P) + wp/2*(Xeff*Xeff);
H = Hq + Hp;
%coherent evolution
spU = -1i* ( spLeftMultiply(H) - spRightMultiply(H) );
%hot bath dissipator
spLh = kh*(nh+1)*(spLrMultiply(kron(sm,idHO)) - 0.5*spLeftMultiply(kron(sp*sm,idHO)) - 0.5*spRightMultiply(kron(sp*sm,idHO))) ...
            + kh*nh*(spLrMultiply(kron(sp,idHO)) - 0.5*spLeftMultiply(kron(sm*sp,idHO)) - 0.5*spRightMultiply(kron(sm*sp,idHO)));

%cold bath dissipator
%use alternative P definition: Approximate by displacement operators
% Idea: p = (1 - exp(i*s*p))/i*s or sin(s*p)/s for s -> 0
% Take sine version for inherent hermitecity. 
% This is +/- displacement by s. Take smallest consistent s, which is dx
Pdx = ( diag(ones(dimp-1,1),1) - diag(ones(dimp-1,1),-1) )/(2*1i*dx);

A = (Xeff+1i*kron(idQ,Pdx))/sqrt(2);
Adg = A';
spLc = kc*(nc+1)*(spLrMultiply(A) - 0.5*spLeftMultiply(Adg*A) - 0.5*spRightMultiply(Adg*A)) ...
            + kc*nc*(spLrMultiply(Adg) - 0.5*spLeftMultiply(A*Adg) - 0.5*spRightMultiply(A*Adg));
%measurement
projP = diag(x>=0)*1.0; %note: 0 with +x
projM = diag(x<0)*1.0;
spLm = km*(spLrMultiply(kron(idQ,projP)) - 0.5*spLeftMultiply(kron(idQ,projP)) - 0.5*spRightMultiply(kron(idQ,projP))) ...
            + km*(spLrMultiply(kron(sx,projM)) - 0.5*spLeftMultiply(kron(idQ,projM)) - 0.5*spRightMultiply(kron(idQ,projM)));
%total liouville op
spLtot = spLh+spLc+spLm+spU;

%% Solve for steady state
sprhoSSAll = spnull(spLtot);
sprhoSS = reshape(sprhoSSAll,dimp*2,dimp*2);
sprhoSS = sprhoSS/trace(sprhoSS);
[rhoSSQ,rhoSSHO] = ptrace(sprhoSS,2,dimp);

%% Get heat flows
Qh = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*H);
Qc = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*H);
W = trace((km*diss(sprhoSS,kron(idQ,projP))+km*diss(sprhoSS,kron(sx,projM)))*H);