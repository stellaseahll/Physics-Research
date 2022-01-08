function [Pp Pm] = getProjector(d)

Pp = diag(ones(1,d)*1/2);
Pm = diag(ones(1,d)*1/2);


for m = 1:d-1
    for n = (m+1):2:d
        M = m-1;
        N = n-1;
        tmp = (-1)^((M+N-1)/2)*2^((M+N)/2)*gamma((M+N)/2+1)*...
            hypergeom([-M,-N],1-M/2-N/2,1/2)/(M+N)/pi/sqrt(factorial(M))/...
            sqrt(factorial(N));
        Pp(m,n) = tmp;
        Pp(n,m) = tmp;
        Pm(m,n) = -tmp;
        Pm(n,m) = -tmp;
    end
end