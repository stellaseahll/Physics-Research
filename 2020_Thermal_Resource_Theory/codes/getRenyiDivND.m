function S = getRenyiDivND(rho,tau,alpha)
%non-diagonal state
[a b] = eig(rho);
idx = find(diag(b)>0)';
lrho = zeros(length(rho));
for i=idx
    lrho = lrho + log(b(i,i))* a(:,i)*a(:,i)';
end
n = length(idx);

[a b] = eig(tau);
idx = find(diag(b)>0)';
ltau = zeros(length(tau));
for i=idx
    ltau = ltau + log(b(i,i))* a(:,i)*a(:,i)';
end

% if nargin = 1, return vn entropy
if (nargin==1)
    S = -trace(rho*(lrho-ltau))/trace(rho);
    return;
end

for i = 1:length(alpha)
    if (alpha(i) == 1)
        S(i) =  trace(rho*(lrho-ltau))/trace(rho);
    elseif (alpha(i) == 0)
        S(i) = n;
%     elseif (alpha(i) < 0)
%         S(i) = log(max(p./q));
    else
        S(i) = log(trace(rho^alpha*tau^(1-alpha))/trace(rho))/(alpha-1);
    end
end