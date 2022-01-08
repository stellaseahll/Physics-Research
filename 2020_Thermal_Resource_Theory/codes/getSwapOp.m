function S = getSwapOp(ds,dp)

S = zeros(ds*dp);

for i = 1:ds
    for j = 1:dp
        S((i-1)*dp+j, (j-1)*dp+i) = 1;
    end
end
