function A = myHermiteCoeff(n)
% Computes matrix of coefficients for normalized Hermite polynomials up to
% order n. Instead of H_n(x) it computes h_n(x) = H_n(x)/sqrt(2^n*n!)
%
%OUTPUT: (n+1)-dim triangular matrix A of coefficients. To evaluate h_k(x) for k=0..n
% use h_k = polyval( A(k+1,1:(k+1)), x )
%INPUT: maximum order n, positive integer
%
%CAUTION: Matrix rows correspond to order of polynomial, columns are
%coefficients of powers of x in DESCENDING order, i.e. A(m+1,1) is
%coefficient of x^m in h_m(x)

A = zeros(n+1,n+1);
A(1,1)=1;
A(2,2)=0;
A(2,1)=sqrt(2);
for nn=2:n % evaluate coefficients c_k for x^k in h_nn
    if nn==2*floor(nn/2) %even nn, all odd k vanish
        A(nn+1,nn+1) = -A(nn,nn-1)/sqrt(2*nn); %k=0
        A(nn+1,1) = sqrt(2/nn)*A(nn,1); %k=nn
        %A(nn+1,2) = sqrt(2/nn)*A(nn,2); %k=nn-1
        for k=2:2:(nn-2)
            A(nn+1,nn+1-k) = sqrt(2/nn)*( A(nn,nn-k+1) - (k+1)/2*A(nn,nn-k-1) );
        end
        
    else %odd nn, all even k vanish
        %A(nn+1,nn+1) = -A(nn,nn-1)/sqrt(2*nn); %k=0
        A(nn+1,1) = sqrt(2/nn)*A(nn,1); %k=nn
        %A(nn+1,2) = sqrt(2/nn)*A(nn,2); %k=nn-1
        for k=1:2:(nn-2)
            A(nn+1,nn+1-k) = sqrt(2/nn)*( A(nn,nn-k+1) - (k+1)/2*A(nn,nn-k-1) );
        end
    end
end
