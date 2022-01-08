clear;clc;
NP = 1; %number of sets measured at one time
nQb = 10000; %total number of qubits
nQbset = 16; %number of qubits per set
nExpt = 50; %average over nExpt experiments
%temperature intervals
Tmin = 1;
Tmax = 10;
nT = 250; %number of intervals
Trange = linspace(Tmin,Tmax,nT+1);
Ttrue = 4;%linspace(Tmin,Tmax,51); %test array of actual temperature 
estT = cell(length(nQbset),length(Ttrue));
for n = 1:length(nQbset) 
    nRun = ceil(nQb/nQbset(n));
    nax{n} = (1:nRun)*nQbset(n);
    d = 2^nQbset(n)-1;
    Eguess = -sqrt(Tmax*Tmin)*log((1/0.85 -1)/d);
%     Eopt = fminsearch(@(E) findOptENEW(E,d,Trange,[],NP),Eguess); %optimal probe for a given Tmax/Tmin
%     pg = TtoProb(Trange,d,Eopt); %probabability of observing ground state for NP 
    estT = cell(1,length(Ttrue));
    for i = 1:length(Ttrue)
        %     estT = zeros(length(E),length(trueT));
        Eguess = -sqrt(Tmax*Tmin)*log((1/0.85 -1)/d);
        Eopt = fminsearch(@(E) findOptENEW(E,d,Trange,[],NP),Eguess); %initial optimal probe for a given Tmax/Tmin independent of prior
        prior = ones(nExpt,1)*findPosterior(Eopt,d,Trange,[],NP,0);
        Eopt = ones(nExpt,1)*Eopt;
        for j = 1:nRun
            fprintf('(%d,%d)\n',i,j);
            pg = 1./(1+d*exp(-Eopt(:,j)*(1./Trange)));
            exppg = 1./(1+d*exp(-Eopt(:,j)./Ttrue(i)));
            expdata = binornd(NP,exppg,nExpt,1);
            likelh = binopdf(expdata*ones(size(Trange)),NP,pg);
            prior = prior.*likelh;
            prior = prior./sum(prior,2);
            estT{n,i}(:,j) = prior*Trange';
            [~,idx] = max(prior');
            est2T{n,i}(:,j) = Trange(idx)';
            meanT{n}(i,j) = mean(estT{i}(:,j));
            for k = 1:nExpt
                Eopt(k,j+1) = fminsearch(@(E) findOptENEW(E,d,Trange,prior(k,:),NP),Eopt(k,j));
            end
        end
        plot(naxis{n},meanT{n}(i,:))
    end
    save('data3.mat');
end
