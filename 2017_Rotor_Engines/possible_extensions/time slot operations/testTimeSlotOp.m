% Test time slot operator (PRE72,066118)
clear;clc;
w = 1;
d = 200;
H = w*(0:d-1);
a = 5;
psi0 = exp(-a^2/2)*(a.^(0:d-1)./sqrt(factorial(0:d-1)));
psi0 = psi0';

t = (0:0.01:1)*2*pi;
t0 = 0;
t1 = pi/4;
psi = zeros(d,length(t));
for i = 1:length(t)
    U = diag(exp(-1i*H*t(i)));
    psi(:,i) = U*psi0;
end

for m = 0:d-1
    for k = 0:d-1
        if (k==m)
            P(k+1,m+1) = exp(-a^2)*a^(2*k)/factorial(k)*(t1-t0);
        else
            P(k+1,m+1) = exp(-a^2)*a^(k+m)/1i/(k-m)/sqrt(factorial(k)*factorial(m))*...
               (exp(-1i*(k-m)*t0)-exp(-1i*(k-m)*t1));
        end
%         if (k==m)
%             P(k+1,m+1) = t1/d;
%         else
%             P(k+1,m+1) = (exp(-1i*w*(k-m)*t1)-1)/(-1i*w*(k-m))/d;
%         end
    end
end
% 
% 
% for m = 0:d-1
%     for k = 0:d-1
%         if ~((abs(k-m)<=m)&&(abs(k-m)<=(d-1)-m))
%             P(k+1,m+1)=0;
%         end
%     end
% end

% for m = 1:d-1
%     P(:,m+1) = P(:,m+1)/norm(P(:,m+1));
%     S = sum((0:d-1)'.*(abs(P(:,m+1)).^2));
%     P(:,m+1) = P(:,m+1)./sqrt(S)*sqrt(m);
% end
% P(:,1) = 0;

for i = 1:length(t)
    U = diag(exp(-1i*H*t(i)));
    psi(:,i) = U*psi0;
    p(i) = sum(abs(P*psi(:,i)).^2);
end
plot(t/2/pi,p)