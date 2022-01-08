function opR = rightMultiply(B,d)
%Converts the Multiplication Y = X*B for (d,d)-matrices into a
%(d^2,d^2) superoperator opL applied like Y(:) = opR*X(:), i.e. to the
%column vector X(:) that vertically stacks the columns of X.
%
%The shape is a matrix of (d,d)-diagonal matrices with the elements of B
%transposed
% [B(1,1)*eye(d), B(2,1)*eye(d), ... B(d,1)*eye(d); ...; B(1,d)*eye(d), ... B(d,d)*eye(d)]
%
%Careful about memory! For 3 qubits: d=8

opR = kron(B.',sparse(eye(d)));