function [ h ] = myHermiteTest (n, x, sw)
%HERMITE Evaulate n-th order NORMALIZED Hermite polynomial H_n(x)/sqrt(2^n*n!)

if sw==0
    h = zeros(1,n+1);
    %n_fact = factorial(n); for unnormalized Hn
    n_fact = sqrt( factorial(n)./2.^n );
    m = 0:floor(n/2);
    h(2*m+1) = n_fact .* (-1).^m ./ (factorial(m) .* factorial(n-2.*m)) .* 2.^(n-2.*m);

    if exist('x','var')
    h = polyval(h, x);
    end 
elseif sw==1
    if (n==0)
        h = 1;
        return;
    end
    if (n==1)
        h = 2*x;
        return;
    end
    Hn(1,1:length(x)) = 1;
    %Hn(2,1:length(x)) = 2*x; %unnormalized Hn
    Hn(2,1:length(x)) = sqrt(2)*x;
    for i = 3:(n+1)
        %Hn(i,1:length(x)) = 2*x.*Hn(i-1,1:length(x))-2*(i-2)*Hn(i-2,1:length(x));
        Hn(i,1:length(x)) = sqrt(2/(i-1))*x.*Hn(i-1,1:length(x)) ...
            - sqrt((i-2)/(i-1))*Hn(i-2,1:length(x));
    end
    h = Hn(end,1:length(x));
elseif sw==2
    A = myHermiteCoeff(n);
    h = polyval(A(n+1,1:(n+1)), x);
end
