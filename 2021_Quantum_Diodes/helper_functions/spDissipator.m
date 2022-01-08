%Computes the usual dissipator D[A](.) = A(.)A' - 1/2* { A'*A, (.) }
% converts full matrices to sparse
function D = spDissipator(A)
    if ~issparse(A)
        A = sparse(A);
    end
    AdA = (A') * A;
    D = spLrMultiply(A) - 1/2*(spLeftMultiply(AdA) + spRightMultiply(AdA));
end