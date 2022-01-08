function M = leftMultiply(A)
    M = (kron(eye(length(A)),A));
end
