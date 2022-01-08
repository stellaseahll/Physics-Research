clear;clc;
NP = 1; %number of sets measured at one time
nQb = 1000; %total number of qubits
nQbset = 1; %number of qubits per set
nExpt = 100; %average over nExpt experiments
%temperature intervals
Tmin = 1;
Tmax = 10;
nT = 200; %number of intervals
Trange = linspace(Tmin,Tmax,nT+1);
Ttrue = 5;%linspace(Tmin,Tmax,51); %test array of actual temperature 
estT = cell(length(nQbset),length(Ttrue));
for n = 1:length(nQbset)

    nRun = ceil(nQb/nQbset(n));
    naxis(n,:) = (1:nRun)*nQbset(n);
    d = 2^nQbset(n)-1;
%     Eopt = fminsearch(@(E) findOptENEW(E,d,Trange,[],NP),Eguess); %optimal probe for a given Tmax/Tmin
%     pg = TtoProb(Trange,d,Eopt); %probabability of observing ground state for NP 
    estT = cell(1,length(Ttrue));
    for i = 1:length(Ttrue)
        %     estT = zeros(length(E),length(trueT));
        xguess = [-sqrt(Tmax*Tmin)*log((1/0.85 -1)) 1];
        xopt = fminsearch(@(x) findOptGround(x,Trange,[],NP),xguess); %initial optimal probe for a given Tmax/Tmin independent of prior
        prior = ones(nExpt,1)*findPosteriorGround(xopt,Trange,[],NP,0);
        xopt = ones(nExpt,1)*xopt;
        gammaopt = xopt(:,2);
        Eopt = xopt(:,1);
        for j = 1:nRun
            fprintf('(%d,%d)\n',i,j);
            nbar = 1./(exp(xopt(:,1)*(1./Trange))-1);
            Gamma = xopt(:,2)*ones(size(Trange)).*(2*nbar+1);
            pg = 1-(nbar.*(1-exp(-Gamma))./(2*nbar+1));
            nbarexp = 1./(exp(xopt(:,1)*(1./Ttrue(i)))-1);
            Gammaexp = xopt(:,2).*(2*nbarexp+1);
            exppg = 1-(nbarexp.*(1-exp(-Gammaexp))./(2*nbarexp+1));
            expdata = binornd(NP,exppg,nExpt,1);
            likelh = binopdf(expdata*ones(size(Trange)),NP,pg);
            prior = prior.*likelh;
            prior = prior./sum(prior,2);
            estT{n,i}(:,j) = prior*Trange';
            [~,idx] = max(prior');
            est2T{n,i}(:,j) = Trange(idx)';
            meanT{n}(i,j) = mean(estT{i}(:,j));
            for k = 1:nExpt
                xopt(k,:) = fminsearch(@(x) findOptGround(x,Trange,prior(k,:),NP),xopt(k,:)); %no save opt values.
            end
            gammaopt(:,j+1) = xopt(:,2); %save opt values
            Eopt(:,j+1) = xopt(:,1);
        end
        plot(naxis(n,:),meanT{n}(i,:))
    end
end