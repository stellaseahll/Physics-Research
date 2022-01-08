function opL = leftMultiply(A,d)
%Converts the Multiplication Y = A*X for (d,d)-matrices into a
%(d^2,d^2) superoperator opL applied like Y(:) = opL*X(:), i.e. to the
%column vector X(:) that vertically stacks the columns of X.
%
%The shape is a matrix with d times A along the diagonal, i.e.
% [A,0...0; 0,A,0...0; ... ; 0...0,A]
%
%Careful about memory! For 3 qubits: d=8
opL = kron(sparse(eye(d)),A);