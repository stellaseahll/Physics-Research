% Class definition for single pass of NP probes interacting via partial
% swap. Thermal bath in between interaction events. Use sparse matrices for definition
% ONLY WORKS IF SYSTEM IS QUBIT
% rho are saved as matrices here, not vectors
%CAUTION:
% Lambda redefined to gamma_t, which is now really gamma*t
% K redefined to G=g*t, which can lead to signs in cosG and sinG
% --> use weird redefined parameters for plots later!
% 
% Redefined SWAP to use XX+YY Hamiltonian!

classdef spSinglePassXYZ < handle
    
    properties
        dS = 2; %system dimension IS FIXED, only implemented for a qubit!
        dP = 2; %probe dimension
        dTot = 4;
        NP = 1; %number of probes
        np = 1; %number of probes in rhoPin
        G = 1; %partial swap parameter, k = sin(gtint) - need to include sign!
        gz = 0; %dephasing parameter
        cosG = cos(1);
        sinG = sin(1);
        gamma_t = 1; %gamma*twait
        nbar = 1; %system temperature
        dnbar = 0.001; %small change for numerical differentiation
        JxS
        JyS
        JzS
        JpS
        JmS
        JxP
        JyP
        JzP
        JpP
        JmP
        JSmat %x,y,z,p,m,mp,pm
        JPmat %x,y,z,p,m,mp,pm
        Umat %all system-probe swap unitaries (matrix form)
%         rhoPin %initial 1-probe state
        rhoPin %input state
        rhoP0 %initial N-probe state
        rhoSP0 %initial state of S+P's
        rhoSS %system steady state
        rhot %final state after NP probes
        rhoAA %final state of first two A
        rhoA %final state of first A
        gsp
        isSteadOp %decides whether system initially in thermal or steady state
        %numerical param for algorithm to compute FI
        % alg: which algorithm (default=1 - mldivide, 2=pcg (needs obj.maxit and obj.tol), 3=lyap)
        alg=2; 
        maxit=1000;
        tol=1e-8;
    end
    
    methods
        function obj = spSinglePassXYZ(g,gz,gammat,nbar,dnbar,rhoP0,np,isSteadOp) 
            obj.nbar = [nbar-dnbar,nbar,nbar+dnbar];
            obj.dnbar = dnbar;
            obj.gamma_t = gammat;
            obj.G = g;
            obj.gz = gz;
            obj.NP = np;
            obj.np = log2(length(rhoP0));
            %check dimension of input state and np
            if (mod(obj.NP,np)~=0) 
                error('Check input state or number of probes');
            end
            obj.rhoPin = rhoP0;
            obj.isSteadOp = isSteadOp;
            obj.dTot = obj.dS*obj.dP^np;
            obj.prepareProbeState();
            obj.getSteadyStates();
            obj.initializeState();
            obj.getOperators(obj.NP);
        end

        function getOperators(obj,np)
            obj.JzS = diag([-0.5,0.5]);
            obj.JpS = [0,0;1,0];
            obj.JmS = obj.JpS';
            obj.JxS = 0.5*(obj.JpS + obj.JmS);
            obj.JyS = -0.5*1i*(obj.JpS - obj.JpS');
            %These are now sparse
            obj.JSmat{1} = kron(obj.JxS,speye(obj.dP^np));
            obj.JSmat{2} = kron(obj.JyS,speye(obj.dP^np));
            obj.JSmat{3} = kron(obj.JzS,speye(obj.dP^np));
            obj.JSmat{4} = kron(obj.JpS,speye(obj.dP^np));
            obj.JSmat{5} = kron(obj.JmS,speye(obj.dP^np));
            obj.JSmat{6} = kron(obj.JmS*obj.JpS,speye(obj.dP^np)); % |g><g|
            obj.JSmat{7} = kron(obj.JpS*obj.JmS,speye(obj.dP^np)); % |e><e|

            JP = (obj.dP-1)/2;
            obj.JzP = diag(-JP:JP);
            m = (-JP:JP).';
            t = sqrt((JP-m).*(JP+m+1));
            obj.JpP = diag(t(1:end-1),-1);
            obj.JmP = obj.JpP';
            obj.JxP = 0.5*(obj.JpP + obj.JmP);
            obj.JyP = -0.5*1i*(obj.JpP - obj.JpP');
            obj.JPmat = cell(5,np);
            JpmP = obj.JpP*obj.JmP; % |e><e|
            JmpP = obj.JmP*obj.JpP; % |g><g|
            %These are also sparse now
            for j = 1:np
                dStoP = obj.dP^(j-1);
                dPtoE = obj.dP^(np-j);
                obj.JPmat{1,j} = kron(speye(obj.dS*dStoP),kron(obj.JxP,speye(dPtoE)));
                obj.JPmat{2,j} = kron(speye(obj.dS*dStoP),kron(obj.JyP,speye(dPtoE)));
                obj.JPmat{3,j} = kron(speye(obj.dS*dStoP),kron(obj.JzP,speye(dPtoE)));
                obj.JPmat{4,j} = kron(speye(obj.dS*dStoP),kron(obj.JpP,speye(dPtoE)));
                obj.JPmat{5,j} = kron(speye(obj.dS*dStoP),kron(obj.JmP,speye(dPtoE)));
                obj.JPmat{6,j} = kron(speye(obj.dS*dStoP),kron(JmpP,speye(dPtoE)));
                obj.JPmat{7,j} = kron(speye(obj.dS*dStoP),kron(JpmP,speye(dPtoE)));
                %get all swap unitaries, make sure it's sparse also when NP=1
                gege = kron( sparse(obj.JmS*obj.JpS), kron( speye(dStoP), kron( JpmP, speye(dPtoE) ) ) );
                egeg = kron( sparse(obj.JpS*obj.JmS), kron( speye(dStoP), kron( JmpP, speye(dPtoE) ) ) );
                egge = kron( sparse(obj.JpS), kron( speye(dStoP), kron( obj.JmP, speye(dPtoE) ) ) );
                %obj.Umat{j} = sin(obj.G)*(egge - egge.'); %U=(G*(|eg><ge| - h.c.))
                obj.Umat{j} = -1i*exp(1i*obj.gz)*sin(obj.G)*(egge + egge.'); % XX+YY swap, U = exp(-1i*G*(|eg><ge| + h.c.))
                obj.Umat{j} = obj.Umat{j} + speye(obj.dS*obj.dP^np); 
                obj.Umat{j} = obj.Umat{j} + (exp(1i*obj.gz)*cos(obj.G)-1)*(gege+egeg);
                if (obj.gz~=0)
                    gggg = kron( sparse(obj.JmS*obj.JpS), kron( speye(dStoP), kron( JmpP, speye(dPtoE) ) ) );
                    eeee = kron( sparse(obj.JpS*obj.JmS), kron( speye(dStoP), kron( JpmP, speye(dPtoE) ) ) );
                    obj.Umat{j} = obj.Umat{j} + (exp(-1i*obj.gz)-1)*(gggg+eeee);
                end
                
            end
        end
        
        function prepareProbeState(obj) %make sparse
            obj.rhoP0=1;
            for i = obj.np:obj.np:obj.NP
                obj.rhoP0 = kron(obj.rhoP0,obj.rhoPin);
            end
        end
        
       
        function initializeState(obj) %should be sparse because obj.rhoP0 is sparse.
            if (~obj.isSteadOp) %start from thermal state of system
                for i = 1:3 %save with slightly shifted nbar for numerical differentiation
                    obj.rhoSP0{i} = kron(obj.prepareThermalState(obj.nbar(i)),obj.rhoP0);
                end
            else %start from steady state of system after many probe interactions
                for i = 1:3
                    obj.rhoSP0{i} = kron(obj.rhoSS{i},obj.rhoP0);
                end
            end
        end
        
        function getSteadyStates(obj) %compute all reduced system steady states for the nbar's
            %This assumes identical qubit ancilla states obj.rhoPin
            if (obj.np ==1) %use the faster speed version
                pA = obj.rhoPin(1,1);
                cA = obj.rhoPin(2,1);
                for i=1:3
                    pth = (obj.nbar(i)+1)/(2*obj.nbar(i)+1);
                    E = exp(-obj.gamma_t*(2*obj.nbar(i)+1));

                    % OBSOLETE: Old swap
                    % A = [ 1-E*cos(obj.G)^2, E*sin(2*obj.G)*abs(cA); ...
                    %       -2*sqrt(E)*sin(obj.G)*abs(cA), 1-sqrt(E)*cos(obj.G) ];
                    % qr = A \ [ E*sin(obj.G)^2*(pA-pth) ; sqrt(E)*sin(obj.G)*abs(cA)*(2*pth-1)];
                    % c = real(qr(2))*exp(1i*phase(cA));

                    % NEW XX+YY swap
                    A = [ 1-E*cos(obj.G)^2, -E*sin(2*obj.G)*abs(cA); ...
                          2*sqrt(E)*sin(obj.G)*abs(cA), 1-sqrt(E)*cos(obj.G) ];
                    qr = A \ [ E*sin(obj.G)^2*(pA-pth) ; -sqrt(E)*sin(obj.G)*abs(cA)*(2*pth-1)];
                    c = 1i*real(qr(2))*exp(1i*phase(cA));

                    p = pth + real(qr(1));
                    obj.rhoSS{i} = [p, conj(c); c, 1-p];
                end
            else
                obj.getOperators(obj.np); %get operators for computing steady state
                options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12,'MaxIterations',1500,'Display','off');
                for i = 1:3    
                    obj.rhoSS{i} = reshape(fsolve(@(rhos) obj.evolution(rhos,i),[1;1;1;1],options),2,2);
                    obj.rhoSS{i} = obj.rhoSS{i}/trace(obj.rhoSS{i});
                end
            end
            
        end
        
        function y = evolution(obj,rhos,idx)
            rho = kron(reshape(rhos,2,2),obj.rhoPin);
            for n = 1:obj.np
                rho = obj.thermalizeSystemEvolve( obj.Umat{n}*(rho*obj.Umat{n}') ,idx);
            end
            rhosf = reshape(obj.ptrace(rho,2,length(obj.rhoPin)),4,1);
            y = rhosf-rhos;          
        end
        
        function rho = thermalizeSystemEvolve(obj,rho0,idx) %thermalizing channel applied to system qubit 
            %total state rho0 -> rho 
            d = obj.dP^obj.np; %total probe dimension

            E = exp(-obj.gamma_t*(2*obj.nbar(idx)+1));
            p = (obj.nbar(idx)+1)/(2*obj.nbar(idx)+1);
            %decay of coherences
            rho = sqrt(E)*rho0;
            %mixing of excited and ground populations
            rho(1:d,1:d) = (p+(1-p)*E) * rho0(1:d,1:d) + p*(1-E) * rho0(d+1:end,d+1:end);
            rho(d+1:end,d+1:end) = (1-p)*(1-E) * rho0(1:d,1:d) + (1-p+p*E) * rho0(d+1:end,d+1:end);
        end
 
        function rho = thermalizeSystem(obj,rho0,idx) %thermalizing channel applied to system qubit 
            %total state rho0 -> rho 
            d = obj.dP^obj.NP; %total probe dimension

            E = exp(-obj.gamma_t*(2*obj.nbar(idx)+1));
            p = (obj.nbar(idx)+1)/(2*obj.nbar(idx)+1);
            %decay of coherences
            rho = sqrt(E)*rho0;
            %mixing of excited and ground populations
            rho(1:d,1:d) = (p+(1-p)*E) * rho0(1:d,1:d) + p*(1-E) * rho0(d+1:end,d+1:end);
            rho(d+1:end,d+1:end) = (1-p)*(1-E) * rho0(1:d,1:d) + (1-p+p*E) * rho0(d+1:end,d+1:end);
        end
 
        
		function rho0 = prepareThermalState(obj,nbar)
            rho0 = diag([nbar+1 nbar]/(2*nbar+1));
        end
		
		function resetState(obj) %just resets the obj.rhot to initial state rhoSP0
            for i=1:3 %init state
                obj.rhot{i} = obj.rhoSP0{i};
            end
		end
		
		function doSwapThermalize(obj,n) %does n-th step of swap and thermalization on all obj.rhot
			for i=1:3
				obj.rhot{i} = obj.thermalizeSystem( obj.Umat{n}*(obj.rhot{i}*obj.Umat{n}') ,i);
			end
		end

		function doIncoherentSwapThermalize(obj,n) %does n-th step of INCOHERENT swap and thermalization on all obj.rhot
			%only makes sense if initial state has no coherences either. Otherwise inconsistent...
			%IF we can always be sure there are NO off-diagonals in S-A states, 
			% then several equivalent defs of Kraus operators for iSWAP are possible.
			% Def 1: 6 Kraus operators, really would kill any coherence { |gg><gg|, |ee><ee|, cosG*|ge><ge|, cosG*|eg><eg|, sinG*|ge><eg|, sinG*|eg><ge| }
			% Def 2: 3 Kraus operators, leaves some coherences intact { |gg><gg| + |ee><ee|, cosG*(|ge><ge| + |eg><eg|), sinG*(|ge><eg| + |eg><ge|) }
			% Both defs differ IF there are coherences, but are consistent for diagonal rho.
			% --> we use Def 2 here, because it's faster. BE SURE YOU KNOW WHAT YOU WANT!
			K1 = obj.JSmat{6}*obj.JPmat{6,n} + obj.JSmat{7}*obj.JPmat{7,n}; %|gg><gg| + |ee><ee|
			K2 = cos(obj.G) * ( obj.JSmat{6}*obj.JPmat{7,n} + obj.JSmat{7}*obj.JPmat{6,n} ); %cosG*(|ge><ge| + |eg><eg|)
			K3 = sin(obj.G) * obj.JSmat{4}*obj.JPmat{5,n};
			K3 = K3+K3';
			for i=1:3
				rhotmp = K1*(obj.rhot{i}*K1) + K2*(obj.rhot{i}*K2) + K3*(obj.rhot{i}*K3);
				obj.rhot{i} = obj.thermalizeSystem( rhotmp ,i);
			end
		end
        
        function getFinalStates(obj) %gets final state for all nbar after NP probes
            %CAUTION: I use ordering: 1st SWAP 2nd Thermalization
            % This is relevant for the initial state, which is the state
            % each probe sees upon interacting
			obj.resetState();
			for n=1:obj.NP
				obj.doSwapThermalize(n);
			end
			%Note that actually, we can remove thermalization from last step,
			%since system will be traced out anyway. But do it for
			%consistency checks. Reduced system state must match steady
			%state then!
        end

        function getIncoherentFinalStates(obj) %gets final state for all nbar w/ incoherent swap
			% Like getFinalStates, but use incoherent swap for Q-Cl comparison. 
			% Executing this function overwrites the final state 
            obj.resetState();
			for n=1:obj.NP
				obj.doIncoherentSwapThermalize(n);
			end
        end

        
        function [F,fz,fy] = getAllFish(obj)
            %computes FI after each step until NP probes, outputs vector of
            %NP FI values (also classical FI, see getFish)
            %To be more efficient and avoid bad scaling, we compute FI only
            %for reduced state of first n probes. So need to ptrace...
            F = [];
			fz = [];
			fy = [];
			obj.resetState();
            for n=1:obj.NP
				obj.doSwapThermalize(n);
                [F(n),fz(n),fy(n)] = obj.getFishSlice(n);
            end
        end
        
        

        function [F,fz,fy] = getNthFish(obj,n)
            %computes FI after each step until NP probes, outputs vector of
            %NP FI values (also classical FI, see getFish)
            %To be more efficient and avoid bad scaling, we compute FI only
            %for reduced state of first n probes. So need to ptrace...
            obj.doSwapThermalize(n);
            [F,fz,fy] = obj.getFishSlice(n);
        end

        function [F,fz,fy] = getAllIncoherentFish(obj)
			%like getAllFish, but using incoherent SWAP.
			F = [];
			fz = [];
			fy = [];
			obj.resetState();
            for n=1:obj.NP
				obj.doIncoherentSwapThermalize(n);
                [F(n),fz(n),fy(n)] = obj.getFishSlice(n);
            end
        end
        
        
        function [F,fz,fy] = getFish(obj)
            %get Q Fisher information F of total state
            %wrapper for getFishSlice (avoid copying code)
            if nargout==3
                [F,fz,fy] = obj.getFishSlice(obj.NP);
            else
                [F,fz] = obj.getFishSlice(obj.NP);
            end
        end
        
        
        
        function [F,fz,fy] = getFishSlice(obj,nP)
            %get Fisher information like in getFish, but only up to nP-th probe
            %need to ptrace then...
            if nP<obj.NP
                rhoP = obj.ptrace( obj.getReducedState(obj.rhot{2}), 2^nP, 2^(obj.NP-nP) );
                drhoP = obj.ptrace( obj.getReducedState(obj.rhot{3} - obj.rhot{1})/2/obj.dnbar, 2^nP, 2^(obj.NP-nP) );
            else
                rhoP = obj.getReducedState(obj.rhot{2});
                drhoP = obj.getReducedState(obj.rhot{3} - obj.rhot{1})/2/obj.dnbar;
            end
            F = obj.catchFish(rhoP,drhoP);
			fz = obj.catchFishClass(diag(rhoP),diag(drhoP));

            %Now for the beam splitter trafo: here we use rotation around y-axis (i.e. in xz plane). 
			% U = sqrt(2)*(Sx + Sz) = [-1,1;1,1]/sqrt(2)
			%CHECK: If we rotate around yz or any rotated vertical plane, do things change?
			%WARNING: Sparsity will be broken, so memory leak for large NP!
			%my tests suggest this kills most FI, not worth it. Maybe wrong choice of UBS...
			if nargout==3 %don't do this shit if not asked for...
				U = 1;
				UBS = sparse([-1,1;1,1]/sqrt(2));
				for n=1:nP %do local trafo on each probe...
					U = kron(U,UBS);
				end
				rhoP = U*(rhoP*U');
				drhoP = U*(drhoP*U');
				fy = obj.catchFishClass(diag(rhoP),diag(drhoP));
			end
        end
        
        function M = getAllSignals(obj,idx)
            %computes triangular NP*NP matrix M with 
            % diagonal entries: Change in n-th probe GS population
            % off-diagonal entries: two-probe correlation in excitation,
            % C=<eg|rho_{ij}|ge>-<e|rho_i|g><g|rho_j|e>
            %Compare magnitudes: Strong quantum enhanced FI should coincide
            %with off-diagonal terms C similar in size as diagonals
            M = zeros(obj.NP);
            for n=1:obj.NP
                %diagonals: population changes
                M(n,n) = trace( (obj.rhot{idx} - obj.rhoSP0{idx}) * obj.JPmat{6,n} ) ;
                for k=1:(n-1) %correlations
                    M(n,k) = trace(obj.rhot{idx}*(obj.JPmat{5,n}*obj.JPmat{4,k})) - ...
                        trace(obj.rhot{idx}*obj.JPmat{5,n})*trace(obj.rhot{idx}*obj.JPmat{4,k});
                end
            end
        end
        
        function rhoP = getReducedState(obj,rho)
            %get reduced state of probes
            %WORKS ONLY IF SYSTEM IS QUBIT!
            d2 = obj.dP^obj.NP;
            rhoP = rho(1:d2,1:d2) + rho(d2+1:end,d2+1:end);
        end
        
        function F = catchFish(obj,rho,drho)
            if obj.alg==1 %do usual FI computation by implicit inverse
                A = kron(rho.',speye(length(rho)))+kron(speye(length(rho)),rho);
                drhot = drho.';
                %Computing inv is very inefficient!
                %F = 2*drhot(:).'*inv(kron(rho.',eye(length(rho)))+kron(eye(length(rho)),rho))*drho(:);
                %mldivide will be very memory inefficient for large matrices!
                %Quite exact, but can go up to 8 or 9 probes only
                F = full( 2*drhot(:).' * ( A \ drho(:) ) );
            elseif obj.alg==2 %use iterative pcg algorithm for F
                %numerical convergence set by parameters maxit and tol, 
                %must be set as object property!
                A = kron(rho.',speye(length(rho)))+kron(speye(length(rho)),rho);
                drhot = drho.';
                %pcg(A,b,tol,maxit) is iterative scheme, will need many iterations at high
                %probe numbers, i.e. long time or else errors in values...
                % may still be more efficient though...
                [x,~,~,~] = pcg(A, drho(:),obj.tol,obj.maxit);
                F = full( 2*drhot(:).' * (x) );
            elseif obj.alg==3 %use lyap method
                %Find FI by solving L*rho+rho*L = 2*drho for L (symmetric),
                % then compute FI = trace(L^2*rho)
                %Don't create a huge superoperator, solve Lyapunov eq directly, 
                %BECAUSE THERE'S A FUCKING MATLAB FUNCTION FOR THIS!!!
                L = lyap(rho,-2*drho);
                F = trace(L*L*rho);
           elseif obj.alg==4 %use iterative lyap method
                %Find FI by solving L*rho+rho*L = 2*drho for L (symmetric),
                % then compute FI = trace(L^2*rho)
                L = lyap2solve(full(rho),full(-2*drho'));
                F = trace(L*L*rho);
            elseif obj.alg==5 %use M.E.S.S. solver
                %B = chol(-2*drho,'lower'); % -2*drho = B*B'
                [Z,D] = mess_lyap(rho,speye(size(rho)),[],-2*drho); %solves A*Z*Z' + Z*Z'*A = -B*B'
                L = Z*D*Z';
                F = trace(L*L*rho);
            else
                error('Count to three cannot? Specify obj.alg=1,2,3!');
            end
        end

		function F = catchFishClass(obj,p,dp)
			%computes classical Fisher info of prob distrib p with numerical derivative dp w.r.t. nbar
			F = sum( (dp.*dp)./p );
		end
        
        function alpha = getExponent(obj,Fn)
            %takes in a vector Fn and assumrs F_n = k (N) ^alpha 
            N = 1:length(Fn);
            c = polyfit(log(N),log(Fn),1);
            alpha = c(1);
        end

        function alpha = getAllExponent(obj,Fn)
            %takes in a vector Fn and assumrs F_n = k N ^alpha 
            for i = 3:length(Fn)
                alpha(i-1) = obj.getExponent(Fn(2:i));
            end
        end

     
        function [A B] = ptrace(obj, AB, dimA, dimB)
            B = zeros(dimB);
            for i = 0:dimA-1
                for j = 0:dimA-1
                    if (i==j)
                        tmp = AB(i*dimB+1:(i+1)*dimB,i*dimB+1:(i+1)*dimB);
                        A(i+1,i+1) = trace(tmp);
                        B = B + tmp;
                    else
                        tmp = AB(i*dimB+1:(i+1)*dimB,j*dimB+1:(j+1)*dimB);
                        A(i+1,j+1) = trace(tmp);
                    end
                end
            end
        end
        
        function F = getProjectorF2(obj,theta)
              Pp = 0.5*eye(2) + (cos(theta)*obj.JzP + sin(theta)*obj.JxP);
              Pm = eye(2)-Pp;
              for i = 1:3
                  pp = trace(obj.rhoAA{i}*kron(Pp,Pp));
                  pm = trace(obj.rhoAA{i}*kron(Pp,Pm));
                  mp = trace(obj.rhoAA{i}*kron(Pm,Pp));
                  mm = trace(obj.rhoAA{i}*kron(Pm,Pm));
                  rho{i} = diag([mm mp pm pp]);
              end
              F = obj.catchFish(rho{2},(rho{3}-rho{1})/obj.dnbar/2);
        end
        
        function [F,theta] = getMaxProjectorF2(obj)
            [theta, F]  = fminsearch(@(theta) -obj.getProjectorF2(theta),1);
            theta = mod(theta,2*pi);
            F = -F;
        end
        
        function getRhoAA(obj)
            obj.rhoAA{1} = obj.ptrace( obj.getReducedState(obj.rhot{1}), 2^2, 2^(obj.NP-2) );
            obj.rhoAA{2} = obj.ptrace( obj.getReducedState(obj.rhot{2}), 2^2, 2^(obj.NP-2) );
            obj.rhoAA{3} = obj.ptrace( obj.getReducedState(obj.rhot{3}), 2^2, 2^(obj.NP-2) );
        end
        
        function getRhoA(obj)
            obj.rhoA{1} = obj.ptrace( obj.getReducedState(obj.rhot{1}), 2^1, 2^(obj.NP-1) );
            obj.rhoA{2} = obj.ptrace( obj.getReducedState(obj.rhot{2}), 2^1, 2^(obj.NP-1) );
            obj.rhoA{3} = obj.ptrace( obj.getReducedState(obj.rhot{3}), 2^1, 2^(obj.NP-1) );
        end
        
        function S = getProjectorDiscord(obj,theta)
%               Pp = 0.5*eye(2) + (cos(theta(1))*obj.JzP + cos(theta(2))*sin(theta(1))*obj.JxP + sin(theta(2))*sin(theta(1))*obj.JyP);
              Pp = 0.5*eye(2) + (cos(theta)*obj.JzP + sin(theta)*obj.JxP);
              Pm = eye(2)-Pp;
              Pp = kron(Pp,eye(2));
              Pm = kron(Pm,eye(2));
              rhop = Pp*obj.rhoAA{2}*Pp;
              probp = trace(rhop);
              rhop = rhop/probp;
              [~,rhobp] = obj.ptrace(rhop,2,2);
              rhom = Pm*obj.rhoAA{2}*Pm;
              probm = trace(rhom);
              rhom = rhom/probm;
              [~,rhobm] = obj.ptrace(rhom,2,2);
              S = probp*obj.getSvn(rhobp) + probm*obj.getSvn(rhobm);
%               [~,rhob] = obj.ptrace(Pp*obj.rhoAA{2}*Pp + Pm*obj.rhoAA{2}*Pm,2,2);
%               S = obj.getSvn(rhob);
        end
        
        function D = getGeometricDiscord(obj)
%               Pp = 0.5*eye(2) + (cos(theta(1))*obj.JzP + cos(theta(2))*sin(theta(1))*obj.JxP + sin(theta(2))*sin(theta(1))*obj.JyP);
            pauli{1} = obj.JxP*2;
            pauli{2} = obj.JyP*2;
            pauli{3} = obj.JzP*2;
            for i = 1:3
                x(i,1) = trace(obj.rhoAA{2}*kron(pauli{i},eye(2)));
                for j = 1:3
                    T(i,j) = trace(obj.rhoAA{2}*kron(pauli{i},pauli{j}));
                end
            end
            T2 = trace(T.'*T);
            x2 = norm(x,2);
            K = x*x.' + T*T.';
            kmax = max(eig(K));
            D = 0.25*(x2+T2-kmax);            
%               [~,rhob] = obj.ptrace(Pp*obj.rhoAA{2}*Pp + Pm*obj.rhoAA{2}*Pm,2,2);
%               S = obj.getSvn(rhob);
        end
        
        function D = getDiscordEnergyBasis(obj)
            SA = obj.getSvn(obj.rhoA{2});
            SAB = obj.getSvn(obj.rhoAA{2});
            S2 = obj.getProjectorDiscord(0);
            D = SA-SAB+S2;
        end
        
        function [D,theta] = getDiscordOpt(obj)
            SA = obj.getSvn(obj.rhoA{2});
            SAB = obj.getSvn(obj.rhoAA{2});
            [theta,S2] = fminsearch(@(theta) obj.getProjectorDiscord(theta),1);
            D = SA-SAB+S2;
        end
        
        
        function S = getSvn(obj,rho)
            k = eig(rho);
            S = -k.*log(k);
            S = sum(S(k>eps));
        end
        
        function C = getConcurrence(obj)
            rhoP = obj.rhoAA{2};
            sy = kron(obj.JyP,obj.JyP);
            rhotmp = sy*rhoP.'*sy;
            R = sqrtm(sqrtm(rhoP)*rhotmp*sqrtm(rhoP));
            K = sort(eig(R),'descend');
            C = max([0, real(K(1)-K(2)-K(3)-K(4))]);
            
        end
        
        function C = getRelCohS(obj,rho)
            %compute relative entropy of coherence
            C = obj.getSvn(diag(diag(rho))) - obj.getSvn(rho); 
        end
        
        function D = getDistributedCoh(obj)
            %arxiv:1801.03919  distributed coherence, Gabriel said its same
            %as discord in energy basis, but let's see.
            D = obj.getRelCohS(obj.rhoAA{2}) - 2*obj.getRelCohS(obj.rhoA{2});
            % 2 becuase reduced state should be the same for both qubits
        end
        
        
        function Fth = getThermalFish(obj)
            Fth = 1/(1+2*obj.nbar(2))^2/(1+obj.nbar(2))/obj.nbar(2);
        end
    end
    
end


