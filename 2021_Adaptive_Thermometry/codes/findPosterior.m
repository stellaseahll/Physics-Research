function post = findPosterior(E,d,Trange,prior,NP,k)
    pg = TtoProb(Trange,d,E); %depends on gap
    pe = 1-pg;
    if isempty(prior)
%         Initial prior = Jeffrey's prior
        tmp = exp(E./Trange);
        prior = d*tmp * E./((d+tmp).^2.*Trange.^2);
        prior(Trange==0) = 0;
        post = prior/sum(prior);
        return;
    end
    likelh = binopdf(k,NP,pg);
    post = likelh.*prior;
    post = post/sum(post);
end