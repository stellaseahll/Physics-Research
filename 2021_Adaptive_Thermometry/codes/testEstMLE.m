clear;clc;
n = 1;
NP = 100;
d = 2^n-1;
Tmin = 1;
Tmax = 10;
nT = 1000; %number of intervals
nExp = 1000;
Trange = linspace(Tmin,Tmax,nT+1);
Tval = 3; %test T values
prior = [];
Eguess = -sqrt(Tmax*Tmin)*log((1/0.85 -1)/d);
Eopt = fminsearch(@(E) findOptMLE(E,d,Trange,prior,NP),Eguess);
pg = TtoProb(Trange,d,Eopt);
estT = cell(1,length(Tval));
for i = 1:length(Tval)
    %     estT = zeros(length(E),length(trueT));
    prior = ones(nExp,1)*findPosterior(Eopt,d,Trange,[],NP,0);
    likelh = binopdf((0:NP)'*ones(size(Trange)),NP,ones(NP+1,1)*pg);
     %number of heads in each experiment

        exppg = 1/(1+d*exp(-Eopt/Tval(i))); %expected ground probability for experiment
        expdata = binornd(NP,exppg,1,nExp);
        prior = prior.*likelh(expdata+1,:);
        prior = prior./sum(prior,2);
        estT = prior*Trange';
%         [~,idx] = max(prior');
%         est2T{i}(:,j) = Trange(idx)';
%         meanT(i) = mean(estT{i}(:,j));
end
