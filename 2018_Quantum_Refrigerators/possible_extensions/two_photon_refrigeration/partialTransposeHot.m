function [mPT,neg,logneg] = partialTransposeHot(m)
%SN 06.02.2018 Computes the partial transpose of a three-qubit matrix
%w.r.t. the first qubit (hot). Use the tensor product convention (ordering)
%[000,001,010,011,100,101,110,111], as follows from kron(H,kron(C,W))
%computes also the negativity and logarithmic negativity if required

%Do manually: (1xy,0ab) <-> (0xy,1ab)
% mPT( 1xy,0ab ) = m( 0xy,1ab ) and mPT( 0xy,1ab ) = m( 1xy,0ab )

mPT = m;

% xy=00, vary ab
mPT(5,1)=m(1,5); mPT(1,5) = m(5,1); % ab=00
mPT(5,2)=m(1,6); mPT(1,6) = m(5,2); % ab=01
mPT(5,3)=m(1,7); mPT(1,7) = m(5,3); % ab=10
mPT(5,4)=m(1,8); mPT(1,8) = m(5,4); % ab=11

% xy=01, vary ab
mPT(6,1)=m(2,5); mPT(2,5) = m(6,1); % ab=00
mPT(6,2)=m(2,6); mPT(2,6) = m(6,2); % ab=01
mPT(6,3)=m(2,7); mPT(2,7) = m(6,3); % ab=10
mPT(6,4)=m(2,8); mPT(2,8) = m(6,4); % ab=11

% xy=10, vary ab
mPT(7,1)=m(3,5); mPT(3,5) = m(7,1); % ab=00
mPT(7,2)=m(3,6); mPT(3,6) = m(7,2); % ab=01
mPT(7,3)=m(3,7); mPT(3,7) = m(7,3); % ab=10
mPT(7,4)=m(3,8); mPT(3,8) = m(7,4); % ab=11

% xy=11, vary ab
mPT(8,1)=m(4,5); mPT(4,5) = m(8,1); % ab=00
mPT(8,2)=m(4,6); mPT(4,6) = m(8,2); % ab=01
mPT(8,3)=m(4,7); mPT(4,7) = m(8,3); % ab=10
mPT(8,4)=m(4,8); mPT(4,8) = m(8,4); % ab=11

if nargout>=2
    %assuming real eigenvals, remove imaginary artefacts
    E = real(eig(mPT));
    neg = sum(abs(E)-E)/2;
    logneg = log2(2*neg+1);
    
end