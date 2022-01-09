function E = getEDiv(rho,beta,H,alpha)
%diagonal state
p = diag(rho);
h = diag(H);
if (nargin==3)
    E = sum(p.*h);
    return;
end

for i = 1:length(alpha)
    if (alpha(i) == 1)
        E(i) =  sum(p.*h);
    else
        E(i) = log(sum(p.^alpha(i).*exp(-beta*h).^(1-alpha(i)))/sum(p.^alpha(i)))/(alpha(i)-1);
    end
end