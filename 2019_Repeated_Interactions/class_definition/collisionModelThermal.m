%Class definition for repeated interaction at Poisson times, system and
%multiple baths. System is arbitrary, bath probes are spins in Gibbs
%states. 

classdef collisionModelThermal < handle
    
    properties
        dS = 2; %system dimension
        NB = 1; %number of baths
        dB = [2]; %array of bath dimensions, length NB 
        TB = [2]; %array of bath temperatures, length NB
        omegaB = [1]; %array of bath gaps, length NB
        gammaB = [1e-2]; %Poisson collision rates per bath
        Gamma %total rate
        whichB %NB array to compute which jump
        dtIntB = [1e-3]; %interaction durations, must be << 1/gamma
        dtSim = 1e-2; %step for simulation
        %NOTE: all of system is stored in eigen basis of HS! baths are
        %given in spin representation (i.e. Jz like HB)
        HS %system Hamiltonian
        Ucomp2E
        %US %trafo matrix of HS eigenbasis
        ES %HS eigenvals as column vect
        ESmES %matrix of -1i*(Ei-Ej)_ij terms for time evolution
        expESmES %matrix of exp(-1i*(Ei-Ej)_ij*dt) terms for time evolution
        rhoS0 %initial system state
        H0 % Free Hamiltonian of bath + system
        Hint %Use baths as left part of tensor product!!!
        Uint %Unitary during interaction, includes free Hamiltonian
        HintEbasis %in total energy eigenbasis HB+HS+Hint
        S % Scattering matrix (defined w.r.t. middle of interaction, unless useVmatrix)
        HB % Bath Hamiltonian
        SB % Entropy of bath 
        UB % Energy of bath 
        USB
        ESB
        ESBmESB
        pthB %thermal Boltzmann weights of bath
        LOps   %all jump operators
        LOpsShort=cell(0); %all jump operators for short time limit
        LOpsScatter %all jump operators using scatterin matrix
        L %total Liouvillian of the master equation valid for long times if free H can be neglected e.g. in 2 qubit scenario
        LScatter %total Liouvillian of the master equation using scattering matrix
        LL % " for one bath
        LLScatter % " for one bath
        LOpseik
        LLeik
        Leik
        Lshort %total Liouvillian of the master equation for short time limit
        rhoSS % steady state
        rhoSSscatter % steady state from scattering matrix
        WSS
        USS
        USSb
        USSs
        SSS
        SSSb
        SSSs
        Seik
        jumpCounter
        tJumpMat
        %by default, use S matrix for master equation. But if manually set
        %to true, then prepareSim() will use V matrix (U0' * U1) instead
        useVmatrix = false
    end
    
    methods
        function obj = collisionModelThermal(ds,Nb,db,ob,Tb,gb,dtb)
            obj.dS = ds;
            obj.NB = Nb;
            obj.HS = zeros(ds);
            %Will take only array of length Nb!
            obj.dB = db(1:Nb);
            obj.TB = Tb(1:Nb);
            obj.omegaB = ob(1:Nb);
            obj.gammaB = gb(1:Nb);
            obj.dtIntB = dtb(1:Nb);
            %system should be initialized in each specific subclass
            %cells of bath operators
            %obj.HB = cell(Nb,1);
            obj.Hint = cell(Nb,1);
            obj.USB = cell(Nb,1);
            obj.ESB = cell(Nb,1);
            obj.ESBmESB = cell(Nb,1);
            %fill variables (OVERRIDE in subclasses)
            obj.rhoS0 = eye(ds)/ds;
            obj.ES = ones(ds,1);
            for n=1:Nb
                obj.Hint{n} = zeros(ds*db(n));
                obj.USB{n} = zeros(ds*db(n));
                obj.ESBmESB{n} = zeros(ds*db(n));
                obj.ESB{n} = zeros(ds*db(n),1);
                pth = exp(-ob(n)/Tb(n) * (0:(db(n)-1)).' );
                pth = pth/sum(pth);
                obj.pthB{n} = pth;
            end
        end
        
        function prepareSim(obj) %preload all stuff for trajectory sim
            %assuming rhoS0 and ES and Hint are set
            %helper matrices for free evolution, precompute for faster loop
            obj.ESmES = -1i*( obj.ES*ones(1,obj.dS) - ones(obj.dS,1)*(obj.ES.') );
            %for interaction
            for n=1:obj.NB
                %NOTE: Bath is left part of tensor product!
                JB = 0.5*( obj.dB(n)-1 );
                obj.HB{n} = diag( (-JB:JB)*obj.omegaB(n) );
                obj.H0{n} = kron(eye(obj.dB(n)),diag(obj.ES)) + kron(obj.HB{n},eye(obj.dS));
                H = obj.Hint{n} + obj.H0{n};

                %S-matrix (or V matrix w.r.t. start point, if manually overriden) !!!
                if obj.useVmatrix
                    obj.S{n} = expm(1i*obj.H0{n}*obj.dtIntB(n)) * expm(-1i*H*obj.dtIntB(n));
                else
                    obj.S{n} = expm(1i*obj.H0{n}*obj.dtIntB(n)/2);
                    obj.S{n} = obj.S{n}*expm(-1i*H*obj.dtIntB(n))*obj.S{n};
                end
                
                [V, h] = eig(H,'vector'); %A*V = V*D
                obj.USB{n} = V;
                obj.HintEbasis{n} = V'*(obj.Hint{n}*V);
                obj.ESB{n} = h;
                %helper matrix for interaction evolution, faster loop (?)
                obj.ESBmESB{n} = -1i*( h*ones(1,obj.dS*obj.dB(n)) - ones(obj.dS*obj.dB(n),1)*(h.') );
				%get Lindblad jump operators, assuming Poisson point process (short interaction times)
                obj.Uint{n} = V * (diag(exp(-1i*obj.dtIntB(n)*h))) *V';
                obj.Seik{n} = expm(-1i*obj.Hint{n}*obj.dtIntB(n));
 				obj.getLOps(n);
                obj.getLOpsEik(n);
                obj.getLOpsScatter(n);
                obj.SB(n) = obj.getEntropy(diag(obj.pthB{n}));
                obj.UB(n) = sum(obj.pthB{n}.*diag(obj.HB{n}));
            end
            %bath drawing: total collision rate, array of cumulative prob
            %distribution for which bath to pick.
            obj.Gamma = sum(obj.gammaB);
            obj.whichB = cumsum(obj.gammaB/obj.Gamma); %probability values between 0 and 1
            %compute Liouvillian of Poisson process and find its steady state
 			obj.getL();
            obj.getLeik();
            obj.getLScatter();
        end
        
        function setdtSim(obj, dt)
            % Call this function to set the timestep for simulation
            %helper matrices for free evolution, precompute for faster loop        
            obj.dtSim = dt;
            obj.expESmES = exp(obj.ESmES*dt);
        end
        
        
        function getLOps(obj,n) 
		%compute Lindblad jump operators for each n-th bath, given the unitary interaction op U.
        %Assumes thermal bath state. 
        %jump op A_{jl} = \sum_ik U_{ijkl} \sqrt{tau_l} |iXk|,U=U_{ijkl} |ijXkl|
		%NOTE: Using this as Lindblad dissipator in ME assumes Poisson point process, i.e. short interaction times!
            obj.LL{n} = zeros(obj.dS^2);
            for j = 1:obj.dB(n)
                for l = 1:obj.dB(n)
                    K = zeros(obj.dB(n));
                    K(j,l) = 1;                     
                    [r,c] = find(kron(K,ones(obj.dS)));
                    tmp = obj.Uint{n}(unique(r),unique(c))*sqrt(obj.pthB{n}(l));
                    obj.LL{n} = obj.LL{n} + obj.lrMultiply(tmp) - ...
                        0.5*( obj.leftMultiply(tmp'*tmp) + obj.rightMultiply(tmp'*tmp));
                    obj.LOps{n,obj.dB(n)*(j-1)+l} = tmp;
                end
            end    
        end
        
       function getLOpsEik(obj,n) 
		%compute Lindblad jump operators for each n-th bath, given the unitary interaction op U.
        %Assumes thermal bath state. 
        %jump op A_{jl} = \sum_ik U_{ijkl} \sqrt{tau_l} |iXk|,U=U_{ijkl} |ijXkl|
		%NOTE: Using this as Lindblad dissipator in ME assumes Poisson point process, i.e. short interaction times!
            obj.LLeik{n} = zeros(obj.dS^2);
            for j = 1:obj.dB(n)
                for l = 1:obj.dB(n)
                    K = zeros(obj.dB(n));
                    K(j,l) = 1;                     
                    [r,c] = find(kron(K,ones(obj.dS)));
                    tmp = obj.Seik{n}(unique(r),unique(c))*sqrt(obj.pthB{n}(l));
                    obj.LLeik{n} = obj.LLeik{n} + obj.lrMultiply(tmp) - ...
                        0.5*( obj.leftMultiply(tmp'*tmp) + obj.rightMultiply(tmp'*tmp));
                    obj.LOpseik{n,obj.dB(n)*(j-1)+l} = tmp;
                end
            end    
        end
        
        function Ldiss = getLdiss(obj) 
		%computes total Liouvillian of a Poisson process, 
		%assuming the bath interaction unitaries as jumps
             Ldiss = zeros(obj.dS^2);
			 %coherent evolution part
             for n = 1:obj.NB
                Ldiss = Ldiss + obj.LLScatter{n}*obj.gammaB(n);
             end
        end
        
        function Hdiss = getHdiss(obj) 
           Hdiss = 1i*( obj.rightMultiply(diag(obj.ES)) - obj.leftMultiply(diag(obj.ES)));
        end
        
        function getLeik(obj)
		%computes total Liouvillian of a Poisson process, 
		%assuming the bath interaction unitaries as jumps
             %obj.L = zeros(obj.dS^2);
			 %coherent evolution part
			 obj.Leik = 1i*( obj.rightMultiply(diag(obj.ES)) - obj.leftMultiply(diag(obj.ES)));
             for n = 1:obj.NB
                obj.Leik = obj.Leik + obj.LLeik{n}*obj.gammaB(n);
             end
        end
        
        function getL(obj) 
		%computes total Liouvillian of a Poisson process, 
		%assuming the bath interaction unitaries as jumps
             %obj.L = zeros(obj.dS^2);
			 %coherent evolution part
			 obj.L = 1i*( obj.rightMultiply(diag(obj.ES)) - obj.leftMultiply(diag(obj.ES)));
             for n = 1:obj.NB
                obj.L = obj.L + obj.LL{n}*obj.gammaB(n);
             end
        end
        function getLOpsScatter(obj,n) 
		%compute Lindblad jump operators for each n-th bath, given the scattering matrix S
        %Assumes thermal bath state. 
        %jump op A_{jl} = \sum_ik U_{ijkl} \sqrt{tau_l} |iXk|,U=U_{ijkl} |ijXkl|
		%NOTE: Using this as Lindblad dissipator in ME assumes Poisson point process, i.e. short interaction times!
            obj.LLScatter{n} = zeros(obj.dS^2);
            for j = 1:obj.dB(n)
                for l = 1:obj.dB(n)
                    K = zeros(obj.dB(n));
                    K(j,l) = 1;                     
                    [r,c] = find(kron(K,ones(obj.dS)));
                    tmp = obj.S{n}(unique(r),unique(c))*sqrt(obj.pthB{n}(l));
                    obj.LLScatter{n} = obj.LLScatter{n} + obj.lrMultiply(tmp) - ...
                        0.5*( obj.leftMultiply(tmp'*tmp) + obj.rightMultiply(tmp'*tmp));
                    obj.LOpsScatter{n,obj.dB(n)*(j-1)+l} = tmp;
                end
            end    
        end
       
        function getLScatter(obj) 
		%computes total Liouvillian of a Poisson process, 
		%assuming the bath interaction unitaries as jumps
             %obj.L = zeros(obj.dS^2);
			 %coherent evolution part
             
			 obj.LScatter = 1i*( obj.rightMultiply(diag(obj.ES)) - obj.leftMultiply(diag(obj.ES)));

             for n = 1:obj.NB
           
                obj.LScatter = obj.LScatter + obj.LLScatter{n}*obj.gammaB(n);
             end
        end
        
        function findSSwithoutSim(obj)
            % run this if have not prepareSim
            obj.prepareSim();
            obj.findSS();
            obj.findSSscatter(); 
        end

        function getSteadyVal(obj)
            for n = 1:obj.NB
                rhoIn = kron(diag(obj.pthB{n}),obj.rhoSS{1});
                rhoOut = obj.Uint{n}*rhoIn* obj.Uint{n}';
                [rhoInB, rhoInS] = ptrace(rhoIn,obj.dB(n),obj.dS);
                [rhoOutB, rhoOutS] = ptrace(rhoOut,obj.dB(n),obj.dS);
                obj.WSS(n) = trace(obj.Hint{n}*(rhoIn-rhoOut));
                obj.USS(n) = trace(obj.H0{n}*(rhoOut-rhoIn));
                obj.USSb(n) = trace(obj.HB{n}*(rhoOutB-rhoInB));
                obj.USSs(n) = trace(diag(obj.ES)*(rhoOutS-rhoInS));
                obj.SSSb(n) = obj.getEntropy(rhoOutB)-obj.getEntropy(rhoInB);
                obj.SSSs(n) = obj.getEntropy(rhoOutS)-obj.getEntropy(rhoInS);
                obj.SSS(n) = obj.getEntropy(rhoOut)-obj.getEntropy(rhoIn);
            end
        end
        
        function S = getEntropy(obj,M)
            v = eig(M);
            S = -sum(v.*log(v));
        end
        
        function findSS(obj) %finds the steady state of the Poisson process Liouvillian L
            %run this if already prepareSim
            
            SS = null(obj.L);
            [~, c] = size(SS);
            if (c==1)
                [v, e] = eig(obj.L,obj.L,'vector','chol');
                id = find(abs(e)<1e-12);
                idx = 1;
                for j = 1:length(id)
                    SS = reshape(v(:,id(j)),[obj.dS,obj.dS]);
                    if (any(diag(SS)<0))
                        continue;
                    end
                    SS = SS/trace(SS);
                    obj.rhoSS{idx} = SS;
                    idx = idx+1;
                end
            else
                for j = 1:c
                    SS = reshape(SS,[obj.dS,obj.dS]);
                    SS = SS/trace(SS);
                    obj.rhoSS{j} = SS;
                end
            end

        end
        
        
        function findSSscatter(obj) %finds the steady state of the Poisson process Liouvillian L
            %run this if already prepareSim
            
            SS = null(obj.LScatter);
            [~, c] = size(SS);
            if (c==0)
   
                [v, e] = eig(obj.LScatter,'vector');
                id = find(abs(e)<1e-12);
                idx = 1;
                for j = 1:length(id)
                    SS = reshape(v(:,id(j)),[obj.dS,obj.dS]);
                    if (any(diag(SS)<0))
                        continue;
                    end
                    SS = SS/trace(SS);
                    obj.rhoSSscatter{idx} = SS;
                    idx = idx+1;
                end
            else
                for j = 1:c
                    SS = reshape(SS,[obj.dS,obj.dS]);
                    SS = SS/trace(SS);
                    obj.rhoSSscatter{j} = SS;
                end
            end

        end
        
        function SSarray = findNull(obj,L) %finds the steady state of any Liouvillian L
            
            SS = null(L);
            [~, c] = size(SS);
            if (c==0)
                [v, e] = eig(L,'vector');
                id = find(abs(e)<1e-12);
                idx = 1;
                for j = 1:length(id)
                    SS = reshape(v(:,id(j)),[obj.dS,obj.dS]);
                    if (any(diag(SS)<0))
                        continue;
                    end
                    SS = SS/trace(SS);
                    obj.rhoSSscatter{idx} = SS;
                    idx = idx+1;
                end
            else
                for j = 1:c
                    SS = reshape(SS,[obj.dS,obj.dS]);
                    SS = SS/trace(SS);
                    SSarray{j} = SS;
                end
            end

        end
        
        function [tWait,bath] = drawJump(obj) %get waiting time for next jump and which bath
            %WARNING: Default is S-matrix formalism, in which case we must
            %draw the MIDDLE POINT of the jump. There we introduce a small
            %cutoff error by assuming jumps cannot overlap (tWait>=tau/2),
            %Here tWait is always the time to START POINT of interaction
            %If V-matrix is used, we draw START POINT directly
            r=rand(1,2);
            tWait = log(1-r(1))/(-obj.Gamma);
            %Find which bath it will be
            bath = sum( r(2)>obj.whichB )+1; %= 1..NB, weighted by relative gamma ratio
            if ~obj.useVmatrix
                tWait = tWait - obj.dtIntB(bath)/2;
                tWait(tWait<0) = eps;
            end
                
            
        end
        
        function [rhoSt, obst] = runMESeik(obj,t,obs)
        % Solves for state at specified time t using Liouvillian L
        % given initial state
            Nt = length(t);
            rhoSt = zeros(obj.dS,obj.dS,Nt);
            rhoSt(:,:,1) = obj.rhoS0;
            rho0 = reshape(obj.rhoS0,obj.dS^2,1);
            Nobs = length(obs);
            obst = zeros(Nobs,Nt); 
            for n=1:Nobs
                obst(n,1) = trace(obs{n}*rhoSt(:,:,1));
            end
            for j = 2:Nt
                rhoSt(:,:,j) = reshape(expm(obj.Leik*t(j))*rho0,obj.dS,obj.dS);
                for n=1:length(obs)
                    obst(n,j) = trace(obs{n}*rhoSt(:,:,j));
                end
            end        
        end

        function [rhoSt, obst] = runMEApprox(obj,t,obs)
        % Solves for state at specified time t using Liouvillian L
        % given initial state
        % Need to get L before running this
            Nt = length(t);
            rhoSt = zeros(obj.dS,obj.dS,Nt);
            rhoSt(:,:,1) = obj.rhoS0;
            rho0 = reshape(obj.rhoS0,obj.dS^2,1);
            Nobs = length(obs);
            obst = zeros(Nobs,Nt); 
            for n=1:Nobs
                obst(n,1) = trace(obs{n}*rhoSt(:,:,1));
            end
            for j = 2:Nt
                rhoSt(:,:,j) = reshape(expm(obj.L*t(j))*rho0,obj.dS,obj.dS);
                for n=1:length(obs)
                    obst(n,j) = trace(obs{n}*rhoSt(:,:,j));
                end
            end        
        end

        function [rhoSt, obst, W] = runMEScatter(obj,t,obs)
        % Solves for state at specified time t using Liouvillian LScatter
        % given initial state
        % HERE: t is array, where t(1)=0 assumed.
        % Need to get LScatter before running this
        % Computes also work powers at each step (which is not a system
        % observable!)
            Nt = length(t);
            rhoSt = zeros(obj.dS,obj.dS,Nt);
            rhoSt(:,:,1) = obj.rhoS0;
            rho0 = reshape(obj.rhoS0,obj.dS^2,1);
            
            %Work powers for each bath
            W = zeros(obj.NB,Nt);
            WOps = cell(obj.NB,1);
            %free unitary needed for S-matrix, but not if V-matrix used
            U0 = cell(obj.NB,1);
            for n=1:obj.NB
                if obj.useVmatrix
                    U0{n} = 1;
                else
                    U0{n} = exp(obj.ESmES*obj.dtIntB(n)/2);
                end
                
                WOps{n} = obj.gammaB(n) * ( obj.S{n}'*( obj.H0{n}*obj.S{n} ) - obj.H0{n} );
                W(n,1) = trace( kron(diag(obj.pthB{n}),U0{n}.*obj.rhoS0)*WOps{n} );
            end
            
            Nobs = length(obs);
            obst = zeros(Nobs,Nt); 
            for n=1:Nobs
                obst(n,1) = trace(obs{n}*rhoSt(:,:,1));
            end
            
            for j = 2:Nt
                rhot = reshape(expm(obj.LScatter*t(j))*rho0,obj.dS,obj.dS);
                rhoSt(:,:,j) = rhot;

                for n=1:obj.NB
                    W(n,j) = trace( kron(diag(obj.pthB{n}),U0{n}.*rhot)*WOps{n} );
                end
                
                for n=1:length(obs)
                    obst(n,j) = trace(obs{n}*rhoSt(:,:,j));
                end
            end
      
        
        end
        
        function [rhoSt, obst, W, Njumps, t] = singleRun(obj,Nt,obs)
            dt = obj.dtSim;
            t = (0:Nt)*dt;
            Njumps = 0;
            %save the whole goddamn state at all goddamn times
            rhoSt = zeros(obj.dS,obj.dS,Nt+1);
            rhoSt(:,:,1) = obj.rhoS0;
            %obs is cell array of observables, given by system operators.
            %Must be given in HS eigenbasis!
            %Will compute trace with system state at all times
            Nobs = length(obs);
            obst = zeros(Nobs,Nt+1); 
            for n=1:Nobs
                obst(n,1) = trace(obs{n}*obj.rhoS0);
            end
            
            %draw waiting time for next jump and decide which bath
%             tJump = -1;
%             while (tJump<0)
%                 [tJump,bath] = obj.drawJump(); 
%                 tJump = tJump-obj.dtIntB(bath)/2; %Random time drawn at center of interaction. Jump starts with displacement of dt/2
%             end
            [tJump,bath] = obj.drawJump(); 
             nJump = floor(tJump/dt)+1; %at which step it happens
            %main loop
            tnow = 0;
            rhoSnow = obj.rhoS0;
            tIntDone = 0;
            W = zeros(1,Nt+1);
            for nt = 1:Nt
                if nt < nJump %no jump occurs in this step
                    rhoSnow = obj.expESmES.*rhoSnow;
                else %jump occurs or is still going on
                    interactionLoop = true;
                    
                    while interactionLoop %necessary, bath interactions may extend over 1 step
                        
                        if tnow < tJump %jump hasn't happened yet
                            %free evolution until jump time
                            rhoSnow = exp(obj.ESmES*(tJump-tnow)).*rhoSnow;
                            tnow=tJump;
							%init operators for bath interaction (only once before jump)
							V = obj.USB{bath}; %basis trafo to eigenbasis of total Hamiltonian
							EmE = obj.ESBmESB{bath}; %gives time evo in eigenbasis representation
							%SB state, but in eigenbasis of total H
							% NOTE: order in kron is B,S ! Makes partial trace easier
                            rhotmp = (V')* (kron(diag(obj.pthB{bath}),rhoSnow) * V);
                            %W(nt) = W(nt)+trace(obj.HintEbasis{bath}*(rhotmp));
                            W(nt) = W(nt)+sum(sum(obj.HintEbasis{bath}.*(rhotmp)));
                        end
						
						%If we arrive here, because of an interaction from the previous step, 
						%count only remaining interaction time. 
                        dtInt = obj.dtIntB(bath)-tIntDone; 
                        
                        if dtInt > nt*dt-tnow %interaction will extend to next step!
                            rhotmp = exp(EmE*(nt*dt-tnow)).*rhotmp; %evolve w/ interaction to end of step
                            %compute reduced system state (must transform back
                            %to product basis and do partial trace)
                            Vtmp = V(1:obj.dS,:);
                            rhoSnow = Vtmp * (rhotmp * (Vtmp') );
                            for k=2:obj.dB(bath)
                                Vtmp = V( ((k-1)*obj.dS+1):(k*obj.dS), : );
                                rhoSnow = rhoSnow + ( Vtmp * (rhotmp * (Vtmp') ) );
                            end
                            %add to interaction time that has been carried out already
							tIntDone = tIntDone + (nt*dt-tnow);
							%now we must proceed to next step, so:
							interactionLoop = false;
                        else %interaction will finish here
                            rhotmp = exp(EmE*dtInt).*rhotmp; %evolve w/ remaining interaction 
                            %compute reduced system state (must transform back
                            %to product basis and do partial trace)
                            %W(nt) = W(nt)-trace(obj.HintEbasis{bath}*(rhotmp));
                            W(nt) = W(nt)-sum(sum(obj.HintEbasis{bath}.*(rhotmp)));
                            Vtmp = V(1:obj.dS,:);
                            rhoSnow = Vtmp * (rhotmp * (Vtmp') );
                            for k=2:obj.dB(bath)
                                Vtmp = V( ((k-1)*obj.dS+1):(k*obj.dS), : );
                                rhoSnow = rhoSnow + ( Vtmp * (rhotmp * (Vtmp') ) );
                            end
                            tnow = tnow + dtInt;
							tIntDone = 0; %reset for next jump
                            Njumps = Njumps + 1; %count the jump 
                            %unfortunately, another jump could occur IN THE SAME FUCKING STEP!
                            [tJump,bath] = obj.drawJump();
                            tJump = tJump + tnow; %time at which next jump happens
                            nJump = floor(tJump/dt)+1; %at which step it happens
                            if nJump>nt %jump will no longer occur in this step
                                %free evolution for remainder
                                rhoSnow = exp(obj.ESmES*(nt*dt-tnow)).*rhoSnow;
                                interactionLoop = false; %break the loop to get to next step
                            end
                        end
                        
                    end
                end
                %At this point, the time step is done, let's save all the
                %shit and get to the next one...
                tnow = nt*dt;
                rhoSt(:,:,nt+1) = rhoSnow; %store the system state
                %store all observable averages the goddamn user wants
                for n=1:Nobs
                    %obst(n,nt+1) = trace(obs{n}*rhoSnow); 
                    obst(n,nt+1) = sum(sum(obs{n}.*rhoSnow)); 
                end
            end
            
        end
        
        function [rhoSt, obst, W, NjumpsAv, t] = manyRuns(obj,Nruns,Nt,obs)
            %runs many trials on single core
            dt = obj.dtSim;
            t = (0:Nt)*dt;
            NjumpsAv = 0;
            rhoSt = zeros(obj.dS,obj.dS,Nt+1);
            %obs is cell array of observables, given by system operators.
            %Must be given in HS eigenbasis!
            %Will compute trace with system state at all times
            Nobs = length(obs);
            obst = zeros(Nobs,Nt+1); 
            W = zeros(1,Nt+1);
            for nr = 1:Nruns
                [rhoSt2, obst2, W2, Njumps] = obj.singleRun(Nt,obs);
                rhoSt = rhoSt + rhoSt2/Nruns;
                obst = obst + obst2/Nruns;
                NjumpsAv = NjumpsAv + Njumps/Nruns;
                W = W + W2/Nruns;
            end
            
        end
        
        function [rhoSt, obst, W, NjumpsAv, t] = manyParRuns(obj,Nruns,Nt,obs)
            %runs many trials on multiple cores (parfor)
            dt = obj.dtSim;
            t = (0:Nt)*dt;
            NjumpsAv = 0;
            rhoSt = zeros(obj.dS,obj.dS,Nt+1);
            %obs is cell array of observables, given by system operators.
            %Must be given in HS eigenbasis!
            %Will compute trace with system state at all times
            Nobs = length(obs);
            obst = zeros(Nobs,Nt+1); 
            W = zeros(1,Nt+1);
            parfor nr = 1:Nruns
                [rhoSt2, obst2, W2, Njumps] = obj.singleRun(Nt,obs);
                rhoSt = rhoSt + rhoSt2/Nruns;
                obst = obst + obst2/Nruns;
                NjumpsAv = NjumpsAv + Njumps/Nruns;
                W = W + W2/Nruns;
            end
            
        end
        
        function [rhoSt, obst, W, S, E, Njumps, t] = singleRunNoise(obj,Nt,obs,sigma)
            dt = obj.dtSim;
            t = (0:Nt)*dt;
            Njumps = 0;
            %save the whole goddamn state at all goddamn times
            rhoSt = zeros(obj.dS,obj.dS,Nt+1);
            rhoSt(:,:,1) = obj.rhoS0;
            %obs is cell array of observables, given by system operators.
            %Must be given in HS eigenbasis!
            %Will compute trace with system state at all times
            Nobs = length(obs);
            obst = zeros(Nobs,Nt+1); 
            for n=1:Nobs
                obst(n,1) = trace(obs{n}*obj.rhoS0);
            end
            
            %draw waiting time for next jump and decide which bath
            [tJump,bath] = obj.drawJump();
            nJump = floor(tJump/dt)+1; %at which step it happens
            %main loop
            tnow = 0;
            rhoSnow = obj.rhoS0;
            tIntDone = 0;
            W = zeros(1,Nt+1);
            S = zeros(3,Nt+1); %row1: entropy of S, 2:bath, 3:total
            E = zeros(3,Nt+1);
            for nt = 1:Nt
                if nt < nJump %no jump occurs in this step
                    rhoSnow = obj.expESmES.*rhoSnow;
                else %jump occurs or is still going on
                    interactionLoop = true;
                    
                    while interactionLoop %necessary, bath interactions may extend over 1 step
                        
                        if tnow < tJump %jump hasn't happened yet
                            %free evolution until jump time
                            rhoSnow = exp(obj.ESmES*(tJump-tnow)).*rhoSnow;
                            S(1,nt) = S(1,nt)-obj.getEntropy(rhoSnow);
                            E(1,nt) = E(1,nt)-sum(obj.ES.*diag(rhoSnow));
                            %frequency is uniform between [(-sigma+1)*omega:
                            %(sigma+1)*omega)
                            ob = (1+(rand()*2*sigma-sigma))*obj.omegaB(bath);
                            pth = exp(-ob/obj.TB(bath) * (0:(obj.dB(bath)-1)).' );
                            pth = pth/sum(pth);
                            obj.pthB{bath} = pth;
                            

                            S(2,nt) = S(2,nt)-obj.getEntropy(diag(obj.pthB{bath}));
                            JB = 0.5*(obj.dB(bath)-1);
                            E(2,nt) = E(2,nt)-ob*sum(pth.*(-JB:JB)');
                            tnow=tJump;
                            %Hamiltonian during interaction between system
                            %and bath, in product basis of bath and energy
                            %eigenbasis of system
                            HBtmp = ob*kron(obj.HB{bath},eye(obj.dS))/obj.omegaB(bath);
                            H0tmp = kron(eye(obj.dB(bath)),diag(obj.ES)) + HBtmp;
                            H = obj.Hint{bath} + H0tmp;
                            rhotmp = kron(diag(obj.pthB{bath}),rhoSnow);       
                            S(3,nt) = S(3,nt) - obj.getEntropy(rhotmp);
                            W(nt) = W(nt)+sum(sum(obj.Hint{bath}.*(rhotmp)));                
                            E(3,nt) = E(3,nt) - sum(sum(H0tmp.*rhotmp));
                        end
						
						%If we arrive here, because of an interaction from the previous step, 
						%count only remaining interaction time. 
                        dtInt = obj.dtIntB(bath)-tIntDone; 
                        
                        if dtInt > nt*dt-tnow %interaction will extend to next step!
                            U = expm(-1i*H*(nt*dt-tnow));
                            rhotmp = U*rhotmp*U'; %evolve w/ interaction to end of step
                            [rhoB, rhoSnow] = ptrace(rhotmp, obj.dB(bath),obj.dS);
                            %add to interaction time that has been carried out already
							tIntDone = tIntDone + (nt*dt-tnow);
							%now we must proceed to next step, so:
							interactionLoop = false;
                        else %interaction will finish here
                            U = expm(-1i*H*(dtInt));
                            rhotmp = U*rhotmp*U'; %evolve w/ remaining interaction 
                            W(nt) = W(nt)-sum(sum(obj.Hint{bath}.*(rhotmp)));
                            E(3,nt) = E(3,nt)+sum(sum(H0tmp.*rhotmp));
                            S(3,nt) = S(3,nt)+obj.getEntropy(rhotmp);
                            [rhoB, rhoSnow] = ptrace(rhotmp, obj.dB(bath),obj.dS);
                            E(2,nt) = E(2,nt) + sum(diag(obj.HB{bath}).*diag(rhoB));
                            E(1,nt) = E(1,nt) + sum((obj.ES).*diag(rhoSnow));
                            S(2,nt) = S(2,nt) + obj.getEntropy(rhoB);
                            S(1,nt) = S(1,nt) + obj.getEntropy(rhoSnow);
                            tnow = tnow + dtInt;
							tIntDone = 0; %reset for next jump
                            Njumps = Njumps + 1; %count the jump 
                            %unfortunately, another jump could occur IN THE SAME FUCKING STEP!
                            [tJump,bath] = obj.drawJump();
                            tJump = tJump + tnow; %time at which next jump happens
                            nJump = floor(tJump/dt)+1; %at which step it happens
                            if nJump>nt %jump will no longer occur in this step
                                %free evolution for remainder
                                rhoSnow = exp(obj.ESmES*(nt*dt-tnow)).*rhoSnow;
                                interactionLoop = false; %break the loop to get to next step
                            end
                        end
                        
                    end
                end
                %At this point, the time step is done, let's save all the
                %shit and get to the next one...
                tnow = nt*dt;
                rhoSt(:,:,nt+1) = rhoSnow; %store the system state
                %store all observable averages the goddamn user wants
                for n=1:Nobs
                    %obst(n,nt+1) = trace(obs{n}*rhoSnow); 
                    obst(n,nt+1) = sum(sum(obs{n}.*rhoSnow)); 
                end
            end
            
        end
        
        function [rhoSt, obst, W, S, E, NjumpsAv, t] = manyRunsNoise(obj,Nruns,Nt,obs,sigma)
            %runs many trials on single core
            dt = obj.dtSim;
            t = (0:Nt)*dt;
            NjumpsAv = 0;
            rhoSt = zeros(obj.dS,obj.dS,Nt+1);
            %obs is cell array of observables, given by system operators.
            %Must be given in HS eigenbasis!
            %Will compute trace with system state at all times
            Nobs = length(obs);
            obst = zeros(Nobs,Nt+1); 
            W = zeros(1,Nt+1);
            E = zeros(3,Nt+1);
            S = zeros(3,Nt+1);
            for nr = 1:Nruns
                [rhoSt2, obst2, W2, S2, E2, Njumps] = obj.singleRunNoise(Nt,obs,sigma);
                rhoSt = rhoSt + rhoSt2/Nruns;
                obst = obst + obst2/Nruns;
                NjumpsAv = NjumpsAv + Njumps/Nruns;
                W = W + W2/Nruns;
                E = E + E2/Nruns;
                S = S + S2/Nruns;
            end
            
        end
        
        function [rhoSt, obst, W, S, E, NjumpsAv, t] = manyParRunsNoise(obj,Nruns,Nt,obs,sigma)
            %runs many trials on multiple cores (parfor)
            dt = obj.dtSim;
            t = (0:Nt)*dt;
            NjumpsAv = 0;
            rhoSt = zeros(obj.dS,obj.dS,Nt+1);
            %obs is cell array of observables, given by system operators.
            %Must be given in HS eigenbasis!
            %Will compute trace with system state at all times
            Nobs = length(obs);
            obst = zeros(Nobs,Nt+1); 
            W = zeros(1,Nt+1);
            E = zeros(3,Nt+1);
            S = zeros(3,Nt+1);
            parfor nr = 1:Nruns
                [rhoSt2, obst2, W2, S2, E2, Njumps] = obj.singleRunNoise(Nt,obs,sigma);
                rhoSt = rhoSt + rhoSt2/Nruns;
                obst = obst + obst2/Nruns;
                NjumpsAv = NjumpsAv + Njumps/Nruns;
                W = W + W2/Nruns;
                E = E + E2/Nruns;
                S = S + S2/Nruns;
            end
            
        end
        
        function M = leftMultiply(obj,A)
            M = (kron(eye(obj.dS),A));
        end
        
        function M = rightMultiply(obj,A)
            M = (kron(A.',eye(obj.dS)));
        end
        
        function M = lrMultiply(obj,A)
            M = (kron(conj(A),A));
        end


    end
    
end