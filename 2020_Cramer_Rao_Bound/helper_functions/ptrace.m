function [A, B] = ptrace(AB, dimA, dimB) %partial trace to get bipartition reduced states

    B = zeros(dimB);
    for i = 0:dimA-1
        for j = 0:dimA-1
            if (i==j)
                tmp = AB(i*dimB+1:(i+1)*dimB,i*dimB+1:(i+1)*dimB);
                A(i+1,i+1) = trace(tmp);
                B = B + tmp;
            else
               tmp = AB(i*dimB+1:(i+1)*dimB,j*dimB+1:(j+1)*dimB);
               A(i+1,j+1) = trace(tmp);
            end
        end 
    end
end