function K = chargeStepsFromEmpty(q,theta,N,Egoal,c)
% Computes the number of charge steps K (charge "time" =K*theta) to
% achieve an ERGOTROPY in the battery >= Egoal, starting from empty. 
% Default = coherent charging
% throws error after at most Kerr steps, so choose Egoal wisely.
% INPUT:
%  q = ground-state population of charge qubits
%  theta = swap angle = g*tau
%  N = battery size (upper most charge level)
%  Egoal = goal energy (in units of gap energy E)
%  c = coherence degree of qubits (DEFAULT = 1). We assume angle alpha=0
%  here.

Kerr = 1000;

if nargin < 5
    c = 1;
end

%define relevant operators and stuff
n = (0:N).'; %column vector of ladder levels (energy in units of gap E)
A = sparse(diag(ones(1,N),1)); %ladder operator (down) sparse
Hd = A + A'; %coherent driving term
P0 = sparse(N+1,N+1); P0(1,1) = 1; %projector on |0>
PN = sparse(N+1,N+1); PN(N+1,N+1) = 1; %projector on |0>

%incoherent terms of charge step trafo
L = sin(theta)^2 * ( q*spDissipator(A) + (1-q)*spDissipator(A') ) ...
    + (1-cos(theta))^2 * ( q*spDissipator(P0) + (1-q)*spDissipator(PN) );

%add coherent terms
L = L - 1i*c*sqrt(q*(1-q))*sin(theta) * ( ...
    cos(theta) * ( spLeftMultiply(Hd) - spRightMultiply(Hd) ) ...
    + (1-cos(theta)) * ( (spLrMultiply2(A,PN)-spLrMultiply2(P0,A)) ...
    + (spLrMultiply2(A',P0)-spLrMultiply2(PN,A')) ) );

rho0 = zeros(N+1,N+1);
rho0(1,1) = 1; %empty initial state
rho = sparse(rho0(:));

E = 0;
%Etmp = 0;
K = 0;

while E<Egoal
    K = K+1;
    %Etmp = E;
    rho = rho + L*rho;
    rhotmp = reshape(rho,N+1,N+1);
    rhotmp = (rhotmp + rhotmp')/2; %remove antihermitian numerical noise part (faster eig)
    E = sum( n .* (diag(rhotmp) - sort( eig( full(rhotmp) ),'descend','ComparisonMethod','real' ) ) ); 
    if K>Kerr
        warning('Egoal cannot be reached!, will stop here')
        break
    end
end

%T = K*theta;



