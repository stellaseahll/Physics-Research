function f = computeIntegral(m,n,a,s)


f = integral(@(x) myHermite(m,x).*myHermite(n,x).*exp(-x.^2).*exp(-(x+a).^2/2/s^2),-inf,inf)/sqrt(2^(m+n)*factorial(m)*factorial(n)*pi); toc;
%
% kend = min([m n]);
% % f = 0;
% % K = factorial(0:max([m,n]));%gamma(1:max([m n])+1);
% H = generateHermite(m+n,a/sqrt(1+2*s^2));
% for k = 0:kend
%     %     A = prod((m-k+1):m)*prod((n-k+1):n)/factorial(k)*2^k/(1+2*s^2)^((m+n)/2-k)
%     A(k+1) = factorial(k)*nchoosek(m,k)*nchoosek(n,k)*2^k/(1+2*s^2)^((m+n)/2-k);
%     B(k+1) = H(m+n-2*k+1);
% %     B(k+1) = hermiteH(m+n-2*k,a/sqrt(1+2*s^2));
%     %     f = f + prod((m-k+1):m)*prod((n-k+1):n)/factorial(k)*2^k/(1+2*s^2)^((m+n)/2-k)*...
%     %         hermiteH(m+n-2*k,a/sqrt(1+2*s^2));
%     %     f = f+A*B;
% end
% f = sum(A.*B);%sum(A(B>0).*B(B>0))+sum(A(B<0).*B(B<0));
% f = f*exp(-a^2/(1+2*s^2))*sqrt((2*s^2)/(1+2*s^2))/(-sqrt(2))^(m+n)/sqrt(gamma(m+1)*gamma(n+1));

% for k = 0:kend
%     f = f + 2^k/K(k+1)/K(m-k+1)/K(n-k+1)/(1+2*s^2)^((m+n)/2-k)*...
%         hermiteH(m+n-2*k,a/sqrt(1+2*s^2));
% end

% f = f*exp(-a^2/(1+2*s^2))*sqrt((2*s^2)/(1+2*s^2))*(-1/sqrt(2))^(m+n)*sqrt(K(m+1)*K(n+1));
% kend = min([m n]);
% f = 0;
% K = gamma(1:(max([m n])+1));
% for k = 0:kend
%     f = f + 2^k/factorial(k)/factorial(m-k)/factorial(n-k)/(1+2*s^2)^((m+n)/2-k)*...
%         hermiteH(m+n-2*k,a/sqrt(1+2*s^2));
% end
%
% f = f*exp(-a^2/(1+2*s^2))*sqrt((2*s^2)/(1+2*s^2))*(-1/sqrt(2))^(m+n)*sqrt(factorial(m)*factorial(n));
% kend = min([m n]);
% f = 0;
% K = factorial(0:max([m n]));%gamma(1:max([m n])+1);
% % K2 = exp(gammaln(m+1)/2+gammaln(n+1)/2-gammaln((0:kend)+1)-gammaln(m-(0:kend)+1)-gammaln(n-(0:kend)+1)+log(2)*((0:kend)-(m+n)/2));
% k = 0:kend;
% K2 = sqrt(K(m+1)*K(n+1))./K(k+1)./K(m-k+1)./K(n-k+1).*(2.^k)./(-sqrt(2))^(m+n)./(1+2*s^2).^((m+n)/2-k)*exp(-a^2/(1+2*s^2))*sqrt((2*s^2)/(1+2*s^2));
% for r = 0:kend
%     f = f + K2(r+1)*hermiteH(m+n-2*r,a/sqrt(1+2*s^2));
% end

end

% 
% 
% function H=generateHermite(maxm,x)
% H(1) = 1;
% H(2) = 2*x;
% for i = 3:(maxm+1)
%     H(i) = 2*x.*H(i-1)-2*(i-2)*H(i-2);
% end
% end