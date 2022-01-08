function errT = findOptGround(x,T,prior,NP)
%x(1) = E, x(2) = gamma
n = 1./(exp(x(1)./T)-1);
Gamma = x(2)*(2*n+1);
pe = n.*(1-exp(-Gamma))./(2*n+1);
pg = 1-pe;
if isempty(prior)
    %  Initial prior = Jeffrey's prior
    num = x(1)^2*(2*exp(x(1)./T).*(-1 + exp(x(2)*coth(x(1)./(2*T))))+(1+exp(x(1)./T)).*x(2).*csch((x(1)./(2*T)).^2)).^2;
    den = 4*(1 + exp(x(1)./T)).^2.*(-1 + exp(x(2)*coth(x(1)./(2*T)))).*(1 + exp(x(1)./T + x(2)*coth(x(1)./(2*T)))).* (T.^4);
    FI = num./den;
    prior = sqrt(FI);
    prior(T==0) = 0;
    prior = prior/sum(prior);
end

likelh = binopdf((0:NP)'*ones(size(T)),NP,ones(NP+1,1)*pg); %NP+1 outcomes (row) * number of true T col
%     size(likelh)
%     size(prior)
post = likelh.*prior;
post = post./sum(post,2);
meanT = post*T';
estT = sum(meanT.*likelh);
errT = sum(abs(T-estT));
end




