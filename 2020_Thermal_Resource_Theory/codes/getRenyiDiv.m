function S = getRenyiDiv(rho,tau,alpha)
%diagonal state
p = diag(rho);
q = diag(tau);
idx = (p>0);
n = sum(p>0);
p = p(idx);
q = q(idx);
% if nargin = 1, return vn entropy
if (nargin==1)
    S = sum(p.*log(p./q));
    return;
end

for i = 1:length(alpha)
    if (alpha(i) == 1)
        S(i) =  sum(p.*log(p./q));
    elseif (alpha(i) == 0)
        S(i) = n;
    elseif (alpha(i) < 0)
        S(i) = log(max(p./q));
    else
        S(i) = log(sum(p.^alpha(i)./(q.^(alpha(i)-1))))/(alpha(i)-1);
    end
end