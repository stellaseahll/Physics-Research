function S = getEntropy(rho,alpha)

p = eig(rho);
p = p(p>0);
% if nargin = 1, return vn entropy
if (nargin==1)
    S = -sum(p.*log(p));
    return;
end

for i = 1:length(alpha)
    if (alpha(i) == 1)
        S(i) =  -sum(p.*log(p));
    elseif (alpha(i) < 0) % alpha-> inf
        S(i) = min(-log(p));
    elseif (alpha(i) == 0)
        S(i) = length(p);
    else
        S(i) = log(sum(p.^alpha(i)))/(1-alpha(i));
    end
end