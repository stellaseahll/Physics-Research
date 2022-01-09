function [P p p2]=plotSwitchHamiltonian(a,t0,t1)

w = 1;
d = 200;
H = w*(0:d-1);
psi0 = exp(-a^2/2)*(a.^(0:d-1)./sqrt(factorial(0:d-1)));
psi0 = psi0';
A = diag(sqrt(1:d-1),1);
t = (0:0.01:1)*2*pi;

psi = zeros(d,length(t));
for i = 1:length(t)
    U = diag(exp(-1i*H*t(i)));
    psi(:,i) = U*psi0;
end

% for m = 0:d-1
%     for k = 0:d-1
%         if (k==m)
%             P(k+1,m+1) = exp(-a^2)*a^(2*k)/factorial(k)*(t1-t0);
%         else
%             P(k+1,m+1) = exp(-a^2)*a^(k+m)/1i/(k-m)/sqrt(factorial(k)*factorial(m))*...
%                (exp(-1i*(k-m)*t0)-exp(-1i*(k-m)*t1));
%         end
% %         if (k==m)
% %             P(k+1,m+1) = t1/d;
% %         else
% %             P(k+1,m+1) = (exp(-1i*w*(k-m)*t1)-1)/(-1i*w*(k-m))/d;
% %         end
%     end
% end

for m = 0:d-1
    for n = 0:d-1
        if (n==(m-1))
            M(n+1,m+1) = exp(-a^2)*a^(2*m)/sqrt(factorial(m)*factorial(m-1))*(t1-t0);
        else
            M(n+1,m+1) = exp(-a^2)*a^(m+n+1)/1i/(m-n-1)/sqrt(factorial(n)*factorial(m))*...
               (exp(-1i*(n-m+1)*t1)-exp(-1i*(n-m+1)*t0));
        end
%         if (k==m)
%             P(k+1,m+1) = t1/d;
%         else
%             P(k+1,m+1) = (exp(-1i*w*(k-m)*t1)-1)/(-1i*w*(k-m))/d;
%         end
    end
end

P = M+M';
P(isnan(P))=0;
for i = 1:length(t)
    U = diag(exp(-1i*H*t(i)));
    psi(:,i) = U*psi0;
    rho = psi(:,i)*psi(:,i)';
    p(i) = trace(P*rho);
    p2(i) = trace((A+A')*rho);
end
p = real(p);
% p = p/max(p);

p2 = real(p2);
% p2 = p2/max(p2);
% clf;
% plot(t/pi,p/max(p));
% hold on;plot(t/pi,p2/max(p2));