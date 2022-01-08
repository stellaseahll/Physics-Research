% function [sprhog,sprhoe,rhoSSQ] = getSteadyEngineXbasisGE(dimp,nh,nc,kh,kc,km,wp,g)
%test measurement engine with x projectors, and HO in xbasis
% clf; 
clear;
Nx = 4;
wp = 1;
ws = 100;
nc = input('nc:');
kc = 0.1; 
kh = 0.001;
km = logspace(-3,-1,50);
g = 2.5;
nh = 1;

% Tc = logspace(-1,0.5,10);
% R = 10;
% dimp = 60;
% wp = 1;
% g = 2.5;
% kc = 0.1; 
% kh = 0.001;
% % km = 0.01;
% ws = 100;
% % Th = 10000; 
% % Tc = 0.1;
% nh = 1; %1/(exp(ws/Th)-1);
% nc = Nc(i);%1/(exp(wp/Tc(i))-1);

%% Operators for qubit
sz = [-1 0; 0 1];
sp = [0 0; 1 0];
sm = sp';
sx = sp+sm;
idQ = eye(2);

%% Operators for HO
%hbar = 1
%mass m = 1
%x = (a+adg)/sqrt(2), so unit should be sqrt(hbar/m/wp)
%p = (a-adg)/1i/sqrt(2)
%x interval from -xmax to xmax
%choice critical! if xmax too high, then wrong results for thermal nbar at
%low dimension dp. If too low, then probably wrong results for high nbar...

dx = g/wp/sqrt(2)/Nx;
xmax = floor(10/dx)*(dx);
x = -xmax:dx:xmax;
dimp = length(x);
p = linspace(-pi/dx,pi/dx,dimp); %fft ordering
N = diag(0:dimp-1);
a = diag(sqrt(1:dimp-1),1);
adg = a';
idHO = eye(dimp);
%position operator in position repres.
X = (diag(x.'));
%FT basis trafo, from position to momentum repres.
FTM = exp(-1i*p.'*x)/sqrt(dimp);
% %check unitarity
% fprintf( 'FT unitarity: %1.2e, %1.2e \n', norm(FTM'*FTM-eye(dimp)),norm(FTM*FTM'-eye(dimp)) );

%momentum operator in position repres.
P = FTM'*(diag(p.')*FTM);
P2 = P*P;
%shifted position operator due to interaction
Xg = X - g/sqrt(2)/wp*idHO;
Xe = X + g/sqrt(2)/wp*idHO;

%% Hamiltonians and Dissipators
% Hq = ws/2*kron(sz,idHO);
Hg = wp/2*P2 + wp/2*(Xg*Xg);
He = wp/2*P2 + wp/2*(Xe*Xe);

%coherent evolution
spUg = -1i* ( spLeftMultiply(Hg) - spRightMultiply(Hg) );
spUe = -1i* ( spLeftMultiply(He) - spRightMultiply(He) );
%hot bath dissipator
spLh1 = kh*(nh+1)*(spLrMultiply(idHO));
spLh2 = kh*nh*(spLrMultiply(idHO));
%cold bath dissipator
Ag = (Xg+1i*P)/sqrt(2);
Agdg = Ag';
AgdgAg = Agdg*Ag;
AgAgdg = Ag*Agdg;
spLcg = kc*(nc+1)*(spLrMultiply(Ag) - 0.5*spLeftMultiply(AgdgAg) - 0.5*spRightMultiply(AgdgAg)) ...
    + kc*nc*(spLrMultiply(Agdg) - 0.5*spLeftMultiply(AgAgdg) - 0.5*spRightMultiply(AgAgdg));
Ae = (Xe+1i*P)/sqrt(2);
Aedg = Ae';
AedgAe = Aedg*Ae;
AeAedg = Ae*Aedg;
spLce = kc*(nc+1)*(spLrMultiply(Ae) - 0.5*spLeftMultiply(AedgAe) - 0.5*spRightMultiply(AedgAe)) ...
    + kc*nc*(spLrMultiply(Aedg) - 0.5*spLeftMultiply(AeAedg) - 0.5*spRightMultiply(AeAedg));
for i = 1:length(km)
    i/length(km)*100
%measurement
projP = diag(x>=0)*1.0; %note: 0 with +x
projM = diag(x<0)*1.0;
spLm1 = km(i)*(spLrMultiply(projM));
spLm2 = spLm1 - km(i)*(spLeftMultiply(projM)+spRightMultiply(projM));

%total liouville op
spLtot = kron([1 0; 0 0],spUg-spLh2+spLcg+spLm2) + kron([0 1; 0 0],spLh1+spLm1)...
   + kron([0 0; 0 1],spUe-spLh1+spLce+spLm2) + kron([0 0; 1 0],spLh2+spLm1);
%% Solve for steady state
sprhoSSAll = spnull(spLtot);
sprhog = reshape(sprhoSSAll(1:dimp^2),dimp,dimp);
sprhoe = reshape(sprhoSSAll((dimp^2+1):(2*dimp^2)),dimp,dimp);
N = trace(sprhog+sprhoe);
% sprhoe = sprhoe/trace(sprhoe);
pg(i) = trace(sprhog)/N;
% rhoSSQ = [trace(sprhog) 0; 0 trace(sprhoe)]/N;
sprhog = sprhog/N;
sprhoe = sprhoe/N;
% plot(x,real(diag(sprhog)),x,real(diag(sprhoe)));
% 
% %% Get heat flows

%% Hamiltonians and Dissipators
% tic;
Hq = ws/2*kron(sz,idHO);
Xeff = kron(idQ,X) + g/sqrt(2)/wp*kron(sz,idHO);
Hp = wp/2*kron(eye(2),P*P) + wp/2*kron(idQ,X*X);
Hint = g*kron(sz,X)/sqrt(2);
sprhoSS = kron([1 0; 0 0],sprhog)+kron([0 0; 0 1],sprhoe);
A = (Xeff+1i*kron(idQ,P))/sqrt(2);
Adg = A';
Qhs(i) = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hq);
Qcs(i) = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hq);
Ws(i) = trace((km(i)*diss(sprhoSS,kron(idQ,projP))+km(i)*diss(sprhoSS,kron(sx,projM)))*Hq);
Qhp(i) = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hp);
Qcp(i) = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hp);
Wp(i) = trace((km(i)*diss(sprhoSS,kron(idQ,projP))+km(i)*diss(sprhoSS,kron(sx,projM)))*Hp);
Qhint(i) = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hint);
Qcint(i) = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hint);
Wint(i) = trace((km(i)*diss(sprhoSS,kron(idQ,projP))+km(i)*diss(sprhoSS,kron(sx,projM)))*Hint);
Qhtot(i) = Qhs(i)+Qhp(i)+Qhint(i);
Qctot(i) = Qcs(i)+Qcp(i)+Qcint(i);
Wtot(i) = Ws(i)+Wp(i)+Wint(i);
end
% X = [0.7:0.1:Xend Xend];
filename = sprintf('Diffgamma_nc%.3f.mat',nc);
save(filename,'nh','nc','kh','kc','g','km','Ws','Qhs','Qcs','Wp','Qhp','Qcp','Wint','Qhint','Qcint','Wtot','Qhtot','Qctot','ws','wp');