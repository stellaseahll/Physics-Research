%Converts left multiplication with matrix A to SPARSE superoperator
function M = spLeftMultiply(A)
    M = kron(speye(size(A)),A); %if one of the matrices is sparse, kron spits out sparse
end
