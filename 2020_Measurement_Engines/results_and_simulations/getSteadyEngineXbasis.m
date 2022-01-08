% function [sprhog,sprhoe,rhoSSQ] = getSteadyEngineXbasis(dimp,nh,nc,kh,kc,km,wp,ws,g)

% function getSteadyEngineXbasis(dimp,nh,nc,kh,kc,km,wp,ws,g)
%test measurement engine with x projectors, and HO in xbasis
clf;
dimp = 50;
wp = 1;
ws = 100;
Th = 100; 
Tc = 0.1; 
nh = 1/(exp(ws/Th)-1);
nc = 1/(exp(wp/Tc)-1);
kh = 0.001;
kc = 0.01;
km = 0.00001;
g = 2;

%% Operators for qubit
sz = [-1 0; 0 1];
sp = [0 0; 1 0];
sm = sp';
sx = sp+sm;
idQ = eye(2);

%% Operators for HO
%hbar = 1
%mass m=1
%x = (a+adg)/spxqrt(2), so unit should be sqrt(hbar/m/wp)
%p = (a-adg)/1i/sqrt(2)
%x interval from -xmax to xmax
%choice critical! if xmax too high, then wrong results for thermal nbar at
%low dimension dp. If too low, then probably wrong results for high nbar...

N = diag(0:dimp-1);
a = diag(sqrt(1:dimp-1),1);
adg = a';
idHO = eye(dimp);
xmax = 12; 
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
A = (Xeff+1i*kron(idQ,P))/sqrt(2);
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
sprhog = sprhoSS(1:dimp,1:dimp);
sprhoe = sprhoSS((1:dimp)+dimp,(1:dimp)+dimp);
rhoSSQ = ptrace(sprhoSS,2,dimp);
px = full([real(sprhog([1 end])) real(sprhoe([1 end]))]);
t1 = FTM*sprhog*FTM';
t2 = FTM*sprhoe*FTM';
pp = full([real(t1([1 end])) real(t2([1 end]))]);
subplot(2,1,1); plot(real(diag(sprhog))); hold on; plot(real(diag(sprhoe))); 
subplot(2,1,2); plot(real(diag(FTM*sprhog*FTM'))); hold on;  plot(real(diag(FTM*sprhoe*FTM'))); 
%% Get heat flows
 tic;
Qh1 = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hq);
Qc1 = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hq);
W1 = trace((km*diss(sprhoSS,kron(idQ,projP))+km*diss(sprhoSS,kron(sx,projM)))*Hq);
Qh2 = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hp);
Qc2 = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hp);
W2 = trace((km*diss(sprhoSS,kron(idQ,projP))+km*diss(sprhoSS,kron(sx,projM)))*Hp);
Qh3 = Qh1+Qh2;
Qc3 = Qc1+Qc2;
W3 = W1+W2;
 toc;