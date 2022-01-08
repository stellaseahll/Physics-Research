clear;clc;
nQb = 10000; %total number of qubits
nQbset = 8; %number of qubits per set
nExpt = 2500; %average over nExpt experiments
NP = 1; %measure NP sets at one go
%temperature intervals
Tmin = 1;
Tmax = 10;
nT = 2500; %number of intervals
Trange = linspace(Tmin,Tmax,nT+1);
Ttrue = linspace(Tmin,Tmax,26); %test array of actual temperature 
estT = cell(length(nQbset),length(Ttrue));
for n = 1:length(nQbset)

    nRun = ceil(nQb/nQbset(n));
    naxis(n,:) = (1:nRun)*nQbset(n);
    d = 2^nQbset(n)-1;
    prior = [];
    Eguess = -sqrt(Tmax*Tmin)*log((1/0.85 -1)/d);
    Eopt = fminsearch(@(E) findOptE(E,d,Trange,prior,NP),Eguess); %optimal probe for a given Tmax/Tmin
    pg = TtoProb(Trange,d,Eopt); %probabability of observing ground state for NP 
    estT = cell(1,length(Ttrue));
    for j = 1:nRun
        %     estT = zeros(length(E),length(trueT));
        prior = ones(nExpt,1)*findPosterior(Eopt,d,Trange,[],NP,0);
        likelh = binopdf((0:NP)'*ones(size(Trange)),NP,ones(NP+1,1)*pg);
        %number of heads in each experiment
        for i = 1:length(Ttrue)
            if mod(j,100)==1
                fprintf('(%d,%d,%d)\n',n,i,j);
            end
            exppg = 1/(1+d*exp(-Eopt/Ttrue(i))); %expected ground probability for experiment
            expdata = binornd(NP,exppg,1,nExpt);
            prior = prior.*likelh(expdata+1,:);
            prior = prior./sum(prior,2);
            estT{n,i}(:,j) = prior*Trange';
            [~,idx] = max(prior');
            est2T{n,i}(:,j) = Trange(idx)';
            meanT{n}(i,j) = mean(estT{i}(:,j));
        end
        save('data.mat');
        plot(naxis(n,:),meanT{n}(i,:)); hold on;
    end
end