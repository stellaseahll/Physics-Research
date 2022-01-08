function M = rightMultiply(A)
    M = (kron(A.',eye(length(A))));
end
