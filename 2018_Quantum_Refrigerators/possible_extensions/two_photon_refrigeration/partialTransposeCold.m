function [mPT,neg,logneg] = partialTransposeCold(m)
%SN 07.02.2018 Computes the partial transpose of a three-qubit matrix
%w.r.t. the 2nd qubit (cold). Use the tensor product convention (ordering)
%[000,001,010,011,100,101,110,111], as follows from kron(H,kron(C,W))
%computes also the negativity and logarithmic negativity if required

%Do manually: (x1y,a0b) <-> (x0y,a1b)
% mPT( x1y,a0b ) = m( x0y,a1b ) and mPT( x0y,a1b ) = m( x1y,a0b )

mPT = m;

% xy=00, vary ab
mPT(3,1)=m(1,3); mPT(1,3) = m(3,1); % ab=00
mPT(3,2)=m(1,4); mPT(1,4) = m(3,2); % ab=01
mPT(3,5)=m(1,7); mPT(1,7) = m(3,5); % ab=10
mPT(3,6)=m(1,8); mPT(1,8) = m(3,6); % ab=11

% xy=01, vary ab
mPT(4,1)=m(2,3); mPT(2,3) = m(4,1); % ab=00
mPT(4,2)=m(2,4); mPT(2,4) = m(4,2); % ab=01
mPT(4,5)=m(2,7); mPT(2,7) = m(4,5); % ab=10
mPT(4,6)=m(2,8); mPT(2,8) = m(4,6); % ab=11

% xy=10, vary ab
mPT(7,1)=m(5,3); mPT(5,3) = m(7,1); % ab=00
mPT(7,2)=m(5,4); mPT(5,4) = m(7,2); % ab=01
mPT(7,5)=m(5,7); mPT(5,7) = m(7,5); % ab=10
mPT(7,6)=m(5,8); mPT(5,8) = m(7,6); % ab=11

% xy=11, vary ab
mPT(8,1)=m(6,3); mPT(6,3) = m(8,1); % ab=00
mPT(8,2)=m(6,4); mPT(6,4) = m(8,2); % ab=01
mPT(8,5)=m(6,7); mPT(6,7) = m(8,5); % ab=10
mPT(8,6)=m(6,8); mPT(6,8) = m(8,6); % ab=11

if nargout>=2
    %assuming real eigenvals, remove imaginary artefacts
    E = real(eig(mPT));
    neg = sum(abs(E)-E)/2;
    logneg = log2(2*neg+1);
    
end