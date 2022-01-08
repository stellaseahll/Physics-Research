function errT = findOptMLE(E,d,trueT,prior,NP)
    % E is an array of 3 numbers
    En = [0 E];
    p = exp(-En'*(1./trueT)); %vector of probabilites depending on temperature
    dT = trueT*0.0001;
    if isempty(prior)
        %         Initial prior = Jeffrey's prior
        pm = exp(-En'*(1./(trueT-dT)));
        pp = exp(-En'*(1./(trueT+dT)));
        dpdT = (pp-pm)/2/dT;
        FI = sum((dpdT).^2./p);
        prior = sqrt(FI);
        prior(trueT==0) = 0;
        prior = prior/sum(prior);
    end

%     likelh = binopdf((0:NP)'*ones(size(trueT)),NP,ones(NP+1,1)*pg); %NP+1 outcomes (row) * number of true T col
%     %     size(likelh)
%     %     size(prior)
%     post = likelh.*prior;
%     post = post./sum(post,2);
%     meanT = post*trueT';
%     estT = sum(meanT.*likelh);
%     errT = sum((log(estT./trueT)).^2);
end