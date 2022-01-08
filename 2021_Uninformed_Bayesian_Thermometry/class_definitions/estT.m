%Biased VS unbiased temperature estimation from simple Gibbs measurements

classdef estTGibbs < handle
    
    properties
        NT %Number of temperature points
        d %Number of degenerate excited states
        optE %Gap between ground and excited states
        T %NT array of temperatures, as given by input
        pT %NT array of excitation prob
        %temperature in units E/kB with E being energy difference between
        %ground and excited
        prior %NT array of prior probabilities for T
        NP %number of qubit probes measured in a run
        %careful: memory! We store NT*(NP+1) matrix of likelihoods
        likh %likelihood function P(n|T) where n is number of excited probes
        lastSample %stores sample from last call of generateData
        lastTrueT %stores true T values of last sample (column vector)
        lastEstT %stores last T estimates (depending which type)
        priorH %function handle for UNNORMALIZED prior as func of p!
        postH %function handle for UNNORMALIZED posterior w/ binomial dist. as func of p!
        postHT %posteriors as func of T (handles)
        normPost %integrated posteriors for all outcomes (over p from 0..1)
    end
    
    methods
        
        function obj = estTGibbs(Tmin,Tmax,nInt,N_P,d,whichPrior)
            %Creates nInt intervals between Tmin and Tmax
            if nargin < 4
                whichPrior = 'Jeffreys'; %by default assume Jeffreys
            end
            obj.d = d;
            obj.T = linspace(Tmin,Tmax,nInt+1)'; %make column vector
            obj.NT = nInt;
            obj.NP = N_P;
            obj.setPrior(whichPrior); 
%             obj.pT = obj.d./(exp(1./obj.T)+obj.d);
%             obj.pT(obj.T==0) = 0;
            
            %likelihood is binomial distribution given excitation prob at
            %NP trials
%             obj.likh = binopdf( ones(obj.NT,1)*(0:obj.NP), obj.NP, obj.pT * ones(1,obj.NP+1) );
        end
        
        function setPrior(obj,whichPrior)
            %sets the prior on the T values specified in array obj.T
            %If numeric value w given, then prior is ~ T^w
            %otherwise can be 'flat' (w=0) or 'Jeffreys'.
            if isnumeric(whichPrior) %assume power law for prior with given power
                obj.prior = obj.T.^whichPrior;
                obj.prior = obj.prior/sum(obj.prior);
                %obj.priorH = @(p) t.^(whichPrior);
            else
                switch whichPrior
                    case 'flat'
                        obj.prior = ones(obj.NT,1)/obj.NT;
                        %obj.priorH = @(t) t.^0;
                    case {'Jeffreys','J','jeffreys','j','jeff','Jeff'}
                        % FI = NP*(dpT/T)^2/(pT*(1-pT)) = NP/T^2 * pT*(1-pT), FI(T=0)=0
                        obj.prior = sqrt(obj.pT.*(1-obj.pT)) ./ obj.T;
                        obj.prior(obj.T==0)=0;
                        obj.prior = obj.prior/sum(obj.prior);
                        %obj.priorH = @(t) sqrt( 1./(exp(1./t)+1).*( 1 - 1./(exp(1./t)+1) ) ) ./ t;
                    otherwise
                        error('Specify valid prior!')
                end
            end
        end
        
        function findOptGap(obj,E)
            obj.pT = obj.d./(exp(E./obj.T)+obj.d);
            obj.pT(obj.T==0) = 0;
            likh = binopdf( ones(obj.NT,1)*(0:obj.NP), obj.NP, obj.pT * ones(1,obj.NP+1) );
            
           	for i = 0:obj.NP
                post(i+1,:) = obj.likh(:,i+1) .* obj.prior;
                post(i+1,:) = post(i+1,:)/sum(post(i+1,:));
                
            end
            
        end
        function generateData(obj,Nruns,trueT)
            %Ttrue can be column vector of length Ntrue, resulting sample
            %will be Ntrue x Nruns array of values from 0..NP
            %if Ttrue is array of different dimension, then it will be made
            %a column vector by Ttrue(:)
            tic
            truep = obj.d./(exp(1./trueT(:))+obj.d);
            obj.lastSample = binornd(obj.NP,truep*ones(1,Nruns));
            obj.lastTrueT = trueT(:);
            toc
        end
        
        function results = pointEstimate(obj)
            %naive point estimate of T based on lastSample (get excitation
            %probability by dividing measured no. of excited probes by NP).
            %Tmean, Tdev give mean, std dev. of T estimates for each Ttrue
            %
            %Careful: This naive point estimator can go out of the
            %specified temperature bounds, so we implement a hard cutoff
            %(which is a bad thing to do, actually, but shouldn't matter
            %too much if not too close to the boundaries)
            tic
            estp = obj.lastSample/obj.NP;
            obj.lastEstT = 1./log(obj.d./estp - obj.d);
            obj.lastEstT(estp==0) = 0;
            
            Tmax = max(obj.T); Tmin = min(obj.T);
            obj.lastEstT(isinf(obj.lastEstT)) = Tmax;
            obj.lastEstT(obj.lastEstT>Tmax) = Tmax;
            obj.lastEstT(obj.lastEstT<Tmin) = Tmin;
            
            results = obj.evalEstimate();
            toc
        end
        
        function results = JesusEstimate(obj)
            %biased point estimator based on 2011.13018
            %must assume equidistant T vals for integration, exclude T=0!
            tic
            [nt,nr] = size(obj.lastSample);
            obj.lastEstT = zeros(nt,nr);
            
            %normalization of likelihood w.r.t. T for estimator
            normL = sum(obj.likh,1); %sum over temperatures
            
            for n=1:nt %loop over trueT values
                %sample contains integer values x=0..NP. With +1, they
                %give column index for likelihood at measured data point x. 
                idx = obj.lastSample(n,:)+1; %array of indices 
                obj.lastEstT(n,:) = exp( log(obj.T).' * ( obj.likh(:,idx)./normL(idx) ) );
            end
            
            results = obj.evalEstimate();
            toc
        end
        
        function results = BayesEstimate(obj)
            %computes Bayes mean estimates from lastSample, individual
            %runs. Memory! We get Nruns*NT array of likelihoods!
            tic
            [nt,nr] = size(obj.lastSample);
            obj.lastEstT = zeros(nt,nr);
            
            for n=1:nt %loop over trueT values
                %sample contains integer values x=0..NP. With +1, they
                %give column index for likelihood at measured data point x. 
                idx = obj.lastSample(n,:)+1; 
                P = obj.likh(:,idx) .* obj.prior;
                P = P ./ sum(P,1);
                obj.lastEstT(n,:) = (obj.T.')*P; %average T from posterior
            end
            
            results = obj.evalEstimate();
            toc
        end
        
        function lp = getPosterior(obj,nruns)
            %computes posterior from lastSample at all Ttrue, but only for
            %the data slice specified by nruns (index 1..Nruns). 
            % In case you want to get an impression of it...
            % --> If nruns is vector of indices, then computes posterior
            %overall posterior by successive Bayesian updates over the
            %slices specified by nruns.
            [nt,~] = size(obj.lastSample);
            lp = obj.prior * ones(1,nt);
            for s=nruns
                idx = obj.lastSample(:,s).'+1;
                lp = lp .* obj.likh(:,idx);
                lp = lp ./ sum(lp,1);
            end
            obj.lastPosterior = lp;         
        end
        
        function results = evalEstimate(obj)
            
            %s = size(obj.lastEstT);
            results = {};
            results{1,1} = 'true T';
            results{1,2} = obj.lastTrueT; 
            results{2,1} = 'mean of estT';
            results{2,2} = mean(obj.lastEstT,2); %average over all Nruns runs
            results{3,1} = 'std of estT';
            results{3,2} = std(obj.lastEstT,0,2); %std dev over all Nruns runs
            results{4,1} = 'mean sqr dev estT 2 trueT';
            results{4,2} = sqrt( mean((obj.lastEstT-obj.lastTrueT).^2,2) );
            %possibly other error cost functions...
        end
        
        
    end
    
end