function U = getSwapU(G,n,N)
%get unitary matrix for partial swap of strength G between 1st and nth
%qubit (out of total N qubits)
%Produces sparse matrix!


% Operators |eg><ge|
gege = sparse([1,0;0,0]);
egeg = sparse([0,0;0,1]);
egge = sparse([0,0;1,0]);
spid = speye(2);


for nP=1:(n-1)
    gege = kron(gege,spid);
    egeg = kron(egeg,spid);
    egge = kron(egge,spid);
end

gege = kron(gege,sparse([0,0;0,1]));
egeg = kron(egeg,sparse([1,0;0,0]));
egge = kron(egge,sparse([0,1;0,0]));

for nP=(n+1):N
    gege = kron(gege,spid);
    egeg = kron(egeg,spid);
    egge = kron(egge,spid);
end

U = sin(G)*(egge - egge.');
U = U + speye(2^(N+1));
U = U + (cos(G)-1)*(gege+egeg);