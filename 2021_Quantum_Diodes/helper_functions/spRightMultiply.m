%Converts right multiplication with matrix A to SPARSE superoperator
function M = spRightMultiply(A)
    M = kron(A.',speye(size(A))); %if one of the matrices is sparse, kron spits out sparse
end
