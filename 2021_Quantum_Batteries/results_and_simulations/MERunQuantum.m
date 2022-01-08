function [cohE, rhof] = MERunQuantum(q,theta,N,tmax,rho0)
%modified from ladderCharging to take input and output

%% Parameters
alpha = 0; %relative phase of qubit superposition (shouldn't matter)
c = 1; %degree of coherence (c=1 superposition, c=0 mixture)

%derived params
p = sin(theta)^2;
v = p*(1-2*q); %drift velocity per step
W = sqrt(q*(1-q))*sin(2*theta); %Rabi parameter
tErg = 2/pi*(p/v^2 - 1); %approx. time from which ergotropy should increase (for incoherent case, when starting at finite charge)
tOpt = N/(c*W+v); %approx. time at which coherent Ergotropy reaches its max (from ground state).
tSat = N/v; %approx. time at which incoherently charging empty battery will saturate at max charge.

%define relevant operators and stuff
n = (0:N).'; %column vector of ladder levels (energy in units of gap E)
A = sparse(diag(ones(1,N),1)); %ladder operator (down) sparse
%Adg = A'; %up
Hd = exp(-1i*alpha)*A; Hd = Hd + Hd'; %coherent driving term
P0 = sparse(N+1,N+1); P0(1,1) = 1; %projector on |0>
PN = sparse(N+1,N+1); PN(N+1,N+1) = 1; %projector on |0>

%sparse superoperators for vectorized Liouville space
%incoherent transformation
Linc = sin(theta)^2 * ( q*spDissipator(A) + (1-q)*spDissipator(A') ) ...
    + (1-cos(theta))^2 * ( q*spDissipator(P0) + (1-q)*spDissipator(PN) );
%coherent trafo
Lcoh = Linc - 1i*sqrt(q*(1-q))*sin(theta) * ( ...
    cos(theta) * ( spLeftMultiply(Hd) - spRightMultiply(Hd) ) ...
    + (1-cos(theta)) * ( exp(-1i*alpha)*(spLrMultiply2(A,PN)-spLrMultiply2(P0,A)) ...
    + exp(1i*alpha)*(spLrMultiply2(A',P0)-spLrMultiply2(PN,A')) ) );

L = c*Lcoh + (1-c)*Linc; %actual trafo at given c

%incoherent random walk: faster if only diagonal elements relevant.
%Propagator matrix for column vector of populations:
K = sin(theta)^2 * ( q*A + (1-q)*A' - speye(N+1));
%must correct boundary terms 
K(1,1) = -sin(theta)^2*(1-q); K(end,end) = -sin(theta)^2*q;
K(1,2) = sin(theta)^2*q; K(end,end-1) = sin(theta)^2*(1-q);

%get passive energy quickly from matrix-form state operator
passE = @(R)(sum( n .* sort( eig( R ),'descend','ComparisonMethod','real' ) ));



%% Check steady states 

% pGibbs = ( (1-q)/q ).^n; pGibbs = pGibbs/sum(pGibbs);

%classical calculation of steady state (DEBUG)

% pSScl = spnull(K); 
% pSScl = pSScl/sum(pSScl);
% ESScl = sum(n.*pSScl);

% rhoSSinc = full(reshape(spnull(Linc),N+1,N+1)); 
%pinc = rhoinc(1:N+2:end); 
% rhoSSinc = rhoSSinc/trace(rhoSSinc);
% pSSinc = diag(rhoSSinc);
%pSSinc = pSSinc/sum(pSSinc);
% ESSinc = sum(n.*pSSinc);
% ErgSSinc = ESSinc - passE(rhoSSinc);

% rhoSScoh = full(reshape(spnull(Lcoh),N+1,N+1));
%pcoh = rhocoh(1:N+2:end); 
% rhoSScoh = rhoSScoh/trace(rhoSScoh);
% pSScoh = diag(rhoSScoh);
%pSScoh = pSScoh/sum(pSScoh);
% ESScoh = sum(n.*pSScoh);
% ErgSScoh = ESScoh - passE(rhoSScoh);
% ErgSSdeph = ESScoh - passE(diag(diag(rhoSScoh)));

%rho = spnull(L);
%p = rho(1:N+2:end); p = p/sum(p);

%% Time evolution coh/incoh/classical

%starting state
% n0 = 0;
% rho0 = zeros(N+1); rho0(n0+1,n0+1) = 1;
p0 = diag(rho0);

% tmax = 5000;
t = (0:tmax);

Ecoh=zeros(1,tmax+1);
Einc=zeros(1,tmax+1);
Ergcoh=zeros(1,tmax+1);
Erginc=zeros(1,tmax+1);
Ergdeph=zeros(1,tmax+1);
Ecl=zeros(1,tmax+1); %for sanity check
trcoh = ones(1,tmax+1);
trinc = ones(1,tmax+1);


%t=0
rcoh = rho0(:);
rinc = rho0(:);
pcl = p0;

Ecoh(1) = sum(n.*rcoh(1:N+2:end));
Einc(1) = sum(n.*rinc(1:N+2:end));
% Ecl(1) = sum(n.*pcl);
Ergcoh(1) = Ecoh(1) - passE(reshape(rcoh,N+1,N+1));
% Erginc(1) = Einc(1) - passE(reshape(rinc,N+1,N+1));
Ergdeph(1) = Ecoh(1) - passE( diag(diag(reshape(rcoh,N+1,N+1))) );
% trcoh(1) = sum(rcoh(1:N+2:end)); %check norm conservation of numerics
% trinc(1) = sum(rinc(1:N+2:end));
% distcoh(1,:) = diag(reshape(rcoh,N+1,N+1)); %distribution 
% distinc(1,:) = diag(reshape(rinc,N+1,N+1)); %distribution 
for k=1:tmax
    
    if mod(k,100)==0
        fprintf('Charge step %i/%i \n',k,tmax);
    end
    
    rcoh = rcoh + Lcoh*rcoh;
    rinc = rinc + Linc*rinc;
    pcl = pcl + K*pcl;
    
    Ecoh(k+1) = sum(n.*rcoh(1:N+2:end));
%     Einc(k+1) = sum(n.*rinc(1:N+2:end));
%     Ecl(k+1) = sum(n.*pcl);
%     distcoh(k+1,:) = diag(reshape(rcoh,N+1,N+1)); %distribution 
%     distinc(k+1,:) = diag(reshape(rinc,N+1,N+1)); %distribution 
    Ergcoh(k+1) = Ecoh(k+1) - passE(reshape(rcoh,N+1,N+1));
%     Erginc(k+1) = Einc(k+1) - passE(reshape(rinc,N+1,N+1));
    Ergdeph(k+1) = Ecoh(k+1) - passE( diag(diag(reshape(rcoh,N+1,N+1))) );
%     trcoh(k+1) = sum(rcoh(1:N+2:end));
%     trinc(k+1) = sum(rinc(1:N+2:end));
    
end

cohE = [Ecoh; Ergcoh; Ergdeph];
% incE = [Einc; Erginc];
rhof = reshape(rcoh,N+1,N+1);

