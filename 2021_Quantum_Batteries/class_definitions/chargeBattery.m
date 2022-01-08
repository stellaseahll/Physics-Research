% Class definition for collisional engine consisting of a qubit as working
% medium and a harmonic oscillator as battery with cutoff dB. The engine
% interacts with thermal probes depending on the control qubits and
% transfers the energy to the harmonic oscillator after
classdef chargeBattery < handle
    
    properties
        dS = 2 %system dimension
        dB = 10 %battery dimension
        dTot = 20
        q = 0.25 %negative temperature used to create system
        battType = 1 %1 = infinite ladder, 2 = finite ladder, 3 = spin
        gt = 0.01 %interaction strength
        %% system operators in S space
        sops
        %% system operators in SB space
        Sops
        %% battery operators in SB space
        Bops
        %% battery operator in B space
        bops
        %% states
        rhoSth %inverted state of system
        rhoB = cell(0) %time evolved reduced battery state
        %% interaction
        V %dimensionless interaction Hamiltonian
        U %unitary
        %% Distribution
        Edist
        %% Energies
        E
        E2
        Erg
    end
    
    methods
        function obj = chargeBattery(q,dS,dB,gt,battType,probeType)
            %% store parameters
            obj.dS = dS;
            obj.dB = dB;
            obj.dTot = dS*dB;
            obj.battType = battType;
            obj.q = q;
            obj.gt = gt;
            %% create states
            obj.createProbeState(probeType);
            obj.rhoB{1} = obj.createGroundBatt(battType);
            obj.Edist(1,:) = diag(obj.rhoB{1});
            %% create operators
            obj.createSops();
            obj.createBops(battType);
            obj.createExchangeInteraction();
        end
        
        function createProbeState(obj,type)
            % create state with negative temp T on a state with dim = dS
%             rat = exp(1/obj.negT); %define rat <1 so we dont get inf
            switch type
                case 'incoh'            
                    obj.rhoSth = obj.createIncoherentProbe(obj.q);
                case 'coh'
                    obj.rhoSth = obj.createCoherentProbe(obj.q);
            end
        end
        
        function rho = createCoherentProbe(obj,q)
            % create state with negative temp T on a state with dim = dS
%             rat = exp(1/obj.negT); %define rat <1 so we dont get inf
            rat = [q, 1-q];
            ket = sqrt(rat)';
            rho = ket*ket';
            rho = rho/trace(rho);
        end
       
        function rho = createIncoherentProbe(obj,q)
            % create state with negative temp T on a state with dim = dS
%             rat = exp(1/obj.negT); %define rat <1 so we dont get inf
            rat = [q, 1-q];
            rho = sparse(diag(rat));
            rho = rho/trace(rho);
        end
        
        function rho = createGroundBatt(obj,battType)
            switch battType
                case 'infladder'
                    rho = zeros(obj.dB);
                    rho(floor(obj.dB/2),floor(obj.dB/2)) = 1;
                    rho = sparse(rho);
                otherwise
                    rho = sparse(diag([1 zeros(1,obj.dB-1)]));
            end
        end
        function createSops(obj)
            % create operators of S in SB space
            % index: 1 = p, 2 = m, 3 = z
            J = (obj.dS-1)/2;
            m = (-J:J).';
            t = sqrt((J-m).*(J+m+1));
            obj.sops{1} = sparse(diag(t(1:end-1),-1));
            obj.sops{2} = obj.sops{1}';
            obj.sops{3} = sparse(diag(m));
            for j = 1:3
                obj.Sops{j} = kron(obj.sops{j},speye(obj.dB));
            end
        end
        
        function createBops(obj,battType)
            % create operators of B in SB space
            % index: 1 = p, 2 = m, 3 = z
            switch battType
                case 'infladder'
                    obj.createBopsInfinLadder()
                case 'finladder'
                    obj.createBopsFinLadder()
                case 'spin'
                    obj.createBopsSpin()
                case 'ho'
                    obj.createBopsHO()
            end
        end
        
        function createBopsInfinLadder(obj)
            m = 0:obj.dB-1;
            obj.bops{1} = sparse(diag(ones(1,obj.dB-1),-1));
            obj.bops{2} = sparse(obj.bops{1}');
            obj.bops{3} = sparse(diag(m));

            for j = 1:3
                obj.Bops{j} = kron(speye(obj.dS),obj.bops{j});
            end
        end
        
        function createBopsFinLadder(obj)
            obj.createBopsInfinLadder(); %DUMMY
%             minN = floor(obj.dB/2)-obj.dB+1;
%             obj.bops{3} = obj.bops{3} - min(diag(obj.bops{3}))*speye(obj.dB);
%             obj.Bops{3} = kron(speye(obj.dS),obj.bops{3});
        end
        
        function createBopsSpin(obj)
            J = (obj.dB-1)/2;
            m = (-J:J).';
            t = sqrt((J-m).*(J+m+1));
            obj.bops{1} = sparse(diag(t(1:end-1),-1));
            obj.bops{2} = sparse(obj.bops{1}');
            obj.bops{3} = sparse(diag(0:obj.dB-1)); %shift energy ground state = 0
            for j = 1:3
                obj.Bops{j} = kron(speye(obj.dS),obj.bops{j});
            end
        end
        
        function createBopsHO(obj)
            obj.bops{3} = sparse(diag(0:obj.dB-1));
            obj.bops{1} = sparse(diag(sqrt(1:obj.dB-1),-1));
            obj.bops{2} = obj.bops{1}';
            for j = 1:3
                obj.Bops{j} = kron(speye(obj.dS),obj.bops{j});
            end
        end
        
        
        function createInteraction(obj,intType)
            %just there in case we want to consider diff interaction
            switch intType
                case 1
                    obj.createExchangeInteraction();
            end
            
        end
        
        function createExchangeInteraction(obj)
            obj.V = (obj.Bops{1}*obj.Sops{2});
            obj.V = obj.V + obj.V';
        end
        
        function runSameInteraction(obj,Nstep)
%             if (nargin == 2)
%                 gt = pi/2;
%             end
%             t = gt/obj.g;
            obj.U = expm(-1i*obj.V*obj.gt);
            for i = 1:Nstep
                obj.runSingleStep();
            end
        end
        
        function runDiffProbes(obj,rho)
            %takes in a cell of rhos that interacts with the battery
            Nstep = length(rho);
            obj.U = expm(-1i*obj.V*obj.gt);
            for i = 1:Nstep
                obj.rhoSth = rho{i};
                obj.runSingleStep();
            end
        end
        
        function runDiffInteraction(obj,gt)
            Nstep = length(gt);
            for i = 1:Nstep
                obj.U = fastExpm(-1i*obj.V*gt(i));
                obj.runSingleStep();
            end
        end
%             
%             
%         end
%         function runFullHamiltonian(obj,Nstep,gt)
%             if (nargin == 2)
%                 gt = pi/2;
%             end
%             t = gt/obj.g;
%             obj.U = expm(-1i*(obj.V+obj.H0)*t);
%             for i = 1:Nstep
%                 obj.runSingleStep();
%             end
%         end
        
        function runSingleStep(obj)
            rhoSB0 = kron(obj.rhoSth,obj.rhoB{end});
            rhoSBf = obj.U*rhoSB0*obj.U';
            [~, obj.rhoB{end+1}] = obj.ptrace(rhoSBf,obj.dS,obj.dB);
            obj.Edist(end+1,:) = diag(obj.rhoB{end});
        end
        
        
        function E = getEnergy(obj)
            for j = 1:length(obj.rhoB)
                E(j) = obj.getExpectation(obj.rhoB{j},obj.bops{3});
            end
        end
        
        
        function E2 = getE2(obj)
            for j = 1:length(obj.rhoB)
                E2(j) = obj.getExpectation(obj.rhoB{j},obj.bops{3}*obj.bops{3});
            end
        end
        
        %% Some math functions
        
        function [A, B] = ptrace(obj, AB, dimA, dimB) %partial trace to get bipartition reduced states
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
            B = sparse(B);
            A = sparse(A);
        end
        
        function E = getExpectation(obj,rho,H)
            %getExpectation computes the expectation value tr(rho*H)
            [a,b] = size(H);
            [c,d] = size(rho);
            if (a~=b)||(c~=d)||(a~=c)
                error('Check input matrices.');
            end
            E = real(full(sum(sum(rho.*(H.'))))); %more efficient than predefined tr(rho*H)
        end
        
        function resetState(obj)
            tmp = obj.rhoB{1};
            clear obj.rhoB
            obj.rhoB = cell(0);
            obj.rhoB{1} = tmp;
        end
        
        function E = getErgotropy(obj)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            H = obj.bops{3};
            for i = 1:length(obj.rhoB)
                rho = full(obj.rhoB{i});
                [pvec, p] = eig(full(rho));
                [p,idx] = sort(diag(p),'descend');
                pvec = pvec(:,idx);
                [evec, e] = eig(full(H));
                E(i) = real(obj.getExpectation(rho,H) - sum(p.*diag(e)));
%                 if (nargout ==2)
%                     U = zeros(length(H));
%                     for k = 1:a
%                         U = U + evec(:,k)*pvec(:,k)';
%                     end
%                     U = sparse(U);
%                 end
            end
        end

        function E = getDephasedErgotropy(obj)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            H = obj.bops{3};
            for i = 1:length(obj.rhoB)
                rho = full(obj.rhoB{i});
                rho = diag(diag(rho));
                [pvec, p] = eig(full(rho));
                [p,idx] = sort(diag(p),'descend');
                pvec = pvec(:,idx);
                [evec, e] = eig(full(H));
                E(i) = real(obj.getExpectation(rho,H) - sum(p.*diag(e)));
%                 if (nargout ==2)
%                     U = zeros(length(H));
%                     for k = 1:a
%                         U = U + evec(:,k)*pvec(:,k)';
%                     end
%                     U = sparse(U);
%                 end
            end
        end
        
        function E = getSystemEnergy(obj)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            H = obj.sops{3};
            E = obj.getExpectation(obj.rhoSth,H);
        end
        
        function S = getSystemEntropy(obj)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            eigVal = eig(obj.rhoSth);
            eigVal = eigVal(eigVal>0);
            S = -sum(eigVal.*log(eigVal));
        end
        
        function E = getSystemErgotropy(obj)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            H = obj.sops{3};
            rho = full(obj.rhoSth);
            [pvec, p] = eig(full(rho));
            [p,idx] = sort(diag(p),'descend');
            pvec = pvec(:,idx);
            [evec, e] = eig(full(H));
            E = obj.getExpectation(rho,H) - sum(p.*diag(e));
%                 if (nargout ==2)
%                     U = zeros(length(H));
%                     for k = 1:a
%                         U = U + evec(:,k)*pvec(:,k)';
%                     end
%                     U = sparse(U);
%                 end
        end
        
        
        function p = getPurity(obj)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            for i = 1:length(obj.rhoB)
                rho =  obj.rhoB{i};
                p(i) = trace(rho*rho);
            end
        end
        
        
    end
    
end

