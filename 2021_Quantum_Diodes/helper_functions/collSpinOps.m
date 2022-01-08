function [jz,jp,jm,jx,jy] = collSpinOps(J)
% Computes sparse matrices for the collective spin operators for given spin
% quantum number J = N/2 (if N is number of qubits combined):
% Jz,J+,J-,Jx,Jy 
% Representation based on Dicke states |J,m> with ascending m=-J..J
% So dimension is 2*J+1

mvals = (-J:1:J).'; %must be column vector
d = 2*J+1;
jpvals = sqrt(J*(J+1) -mvals.*(mvals+1) );

jz = sparse(d,d);
jp = sparse(d,d);

jz = spdiags(mvals,0,jz); %places mvals into main diagonal of sparse mat

jp = spdiags(jpvals(1:end-1),-1,jp); %1st lower co-diagonal

jm = jp';
jx = (jp+jm)/2;
jy = -1i/2*(jp-jm);

