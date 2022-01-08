function [W,Qh,Qc,sprhoSS] = getSteadyEngineXbasis_POVM(dimp,nh,nc,kh,kc,km,wp,ws,g,z0,sigma)
%test measurement engine with Gaussian POVM (of width sigma), where bitflip
%is applied if result z < z0, otherwise no bitflip. 
%This should converge to projective x-measurement for sigma -> 0
%But might be useful when using with large sigma, because we expect better 
%numerical convergence behavior due to smoothness in position.
%
%superoperator can be constructed explicitly in position representation,
%and is quite sparse because diagonal (apart from bitflip)
%
%Standard choice: take sigma of oscillator GS for a reasonable precision,
% take z0=0 to set threshold at origin.
% But engine might perform better for z0 < 0, at least measurement should
% be more selective then.
%Sanity check: Should get abt same result as projective case for z0=0 and
% sigma -> 0 (i.e. << position resolution)


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

%% create superoperators for POVM
%diagonal elements of operator kron(idQ,x)
%tic
id_x = kron([1;1],x.'); %vector in sz-x representation

%diagonal elements of left or right multiplication with id_x
%keep it dense, because we need to apply functions later also on zero-elements
XL = kron(ones(2*dimp,1),id_x); %diagonal elements of leftMultiply of id_x
XR = kron(id_x,ones(2*dimp,1)); %diagonal elements of rightMultiply of id_x

%position-dependent functions f(x1,x2) to multiply to <x1|rho|x2>
%these superoperators are diagonal in our sz-x representation, so compute only
%the vectors of diagonal elements first
gM = exp(-(XL-XR).^2/8/sigma^2);
erfP = erf( (z0-(XL+XR)/2)/(sqrt(2)*sigma) );
f_flip = gM.*(1+erfP)/2; %function to apply together with bitflip sx * rho *sx
f_noflip = gM.*(1-erfP)/2 - 1; %to apply when no flip, contains the - rho part

%Now make them actual superoperators. For this need to create diagonal
%sparse matrix from diagonal elements
iv = (1:length(gM)).'; %need vector of indices in superoperator space
F_flip = sparse( iv, iv, f_flip ); % F_flip(iv(k),iv(k))=f_flip(k)
F_noflip = sparse( iv, iv, f_noflip );
bitflip = spLrMultiply(kron(sx,idHO)); %just the bitflip operation

%total measurement superoperator in sz-x repres
% the matrix multiplication can probably be done faster somehow...
spLm = km * ( F_noflip + bitflip*F_flip );
%actual formula is:
% <x1|Lm*rho|x2> = km*{ gM(x1,x2)*[ (1+erfP(x1,x2))/2 * <x1|sx*rho*sx|x2> + 
%                     (1-erfP(x1,x2))/2 * <x1|rho|x2> ] - <x1|rho|x2> }
%
%Even if I messed up the ordering of left/right multiply, shouldn't matter,
%as all functions are symmetric under x1 <--> x2
%whos spLm
%toc

%% Hamiltonians and Dissipators
%tic
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

%total liouville op
spLtot = spLh+spLc+spLm+spU;
%whos spLtot
%toc

%% Solve for steady state
%tic
sprhoSSAll = spnull(spLtot);
sprhoSS = reshape(sprhoSSAll,dimp*2,dimp*2);
sprhoSS = sprhoSS/trace(sprhoSS);
[rhoSSQ,rhoSSHO] = ptrace(sprhoSS,2,dimp);
%toc

%% Get heat flows
Qh = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*H);
Qc = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*H);

LmRhoSS = reshape(spLm*sprhoSS(:),dimp*2,dimp*2);
W = trace(LmRhoSS*H);