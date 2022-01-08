%Converts left-right multiplication with SPARSE matrix A, A*rho*A', to SPARSE superoperator
function M = spLrMultiply(A)
    if ~issparse(A)
        A = sparse(A);
    end
    M = (kron(conj(A),A));
end
