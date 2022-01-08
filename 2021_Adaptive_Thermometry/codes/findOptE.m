function errT = findOptE(E,d,trueT,prior,NP)
    pg = TtoProb(trueT,d,E); %depends on gap
    pe = 1-pg;
    if isempty(prior)
%         Initial prior = Jeffrey's prior

        prior = sqrt(pg.*pe) ./ trueT;
        prior(trueT==0) = 0;
        prior = prior/sum(prior);
    end

    likelh = binopdf((0:NP)'*ones(size(trueT)),NP,ones(NP+1,1)*pg); %NP+1 outcomes (row) * number of true T col
%     size(likelh)
%     size(prior)
    post = likelh.*prior;
    post = post./sum(post,2);
    meanT = post*trueT'
    
    estT = zeros(size(trueT));
     likelh(:,1)
%     prior = ones(NP+1,1)*prior;
    for k = 1:length(trueT)
        estT(k) = estT(k) + meanT'*likelh(:,k);
    end
    errT = sum(abs(trueT-estT));
end