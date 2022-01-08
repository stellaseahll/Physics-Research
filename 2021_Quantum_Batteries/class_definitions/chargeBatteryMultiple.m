% Class definition for collisional engine consisting of a qubit as working
% medium and a harmonic oscillator as battery with cutoff dB. The engine
% interacts with thermal probes depending on the control qubits and
% transfers the energy to the harmonic oscillator after
classdef chargeBatteryMultiple < handle
    
    properties
        dS = 2 %system dimension
        nS = 1 %default 1 battery
        dB = 10 %battery dimension
        dTot = 20
        negT = -1 %negative temperature used to create system
        battType = 0 %0 = infinite ladder, 1 = finite ladder, 2 = spin
        g = 0.01 %interaction strength
        %% system operators in S space
        sops
        %% system operators in SB space
        Sops
        %% battery operators in SB space
        Bops
        %% battery operator in B space
        bops
        %% states
        rhoIndSth
        rhoSth %inverted state of system
        rhoSf = cell(0)
        rhoB = cell(0) %time evolved reduced battery state
        rhoBS0 %initial state
        %% energies
        Eth
        Ergth
        EB
        ErgB
        
        %% interaction
        HS0 %bare Hamiltonian
        parV %parallel interaction Hamiltonian
        twoBodyV = cell(0) %two body interaction between probes and battery
        U %unitary
        
    end
    
    methods
        function obj = chargeBatteryMultiple(negT,dS,dB,nS,battType)
%             if (negT>0)
%                 error('Negative temperature must be negative.');
%             end
            
            %% store parameters
            obj.dS = dS;
            obj.dB = dB;
            obj.nS = nS;
            obj.dTot = dS*dB;
            obj.battType = battType;
            obj.negT = negT;
            %% create operators
            obj.createOperators;
            %% create states
            obj.createInvTState();
            obj.createGroundBatt();
            %             obj.createInitState();
%             %% energies
%             obj.getInvTEnergy();
        end
        
        function rho = createCoherentState(obj)
            prob = (exp(1/obj.negT)).^(obj.dS:-1:1);
            prob = prob/sum(prob);
            rho = sqrt(prob)'*sqrt(prob);
        end
        
        
        
        function rho = createInvTState(obj)
            % create state with negative temp T on a state with dim = dS
            rat = exp(1/obj.negT); %define rat <1 so we dont get inf.
            rho = sparse(diag(rat.^(obj.dS:-1:1)));
            rho = rho/trace(rho);
        end
        
        function rho = createGroundBatt(obj)
            if (obj.battType == 2 || obj.battType ==3)
                rho = sparse(diag([1 zeros(1,obj.dB-1)]));
            elseif (obj.battType == 1)
                rho = zeros(obj.dB);
                rho(ceil(obj.dB/4),ceil(obj.dB/4)) = 1;
                rho = sparse(rho);
            end
            obj.rhoB{1} = rho;
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
            %             for j = 1:3
            %                 obj.Sops{j} = kron(obj.sops{j},speye(obj.dB));
            %             end
        end
        
        function createBops(obj)
            % create operators of B in SB space
            % index: 1 = p, 2 = m, 3 = z
            switch obj.battType
                case 1
                    obj.createBopsInfinLadder()
                case 2
                    obj.createBopsFinLadder()
                case 3
                    obj.createBopsSpin()
            end
        end
        
        function createBopsInfinLadder(obj)
            maxN = floor(obj.dB/4*3);
            m = maxN-obj.dB+1:maxN;
            obj.bops{1} = sparse(diag(ones(1,obj.dB-1),-1));
            obj.bops{2} = sparse(obj.bops{1}');
            obj.bops{3} = sparse(diag(m));
            %             for j = 1:3
            %                 obj.Bops{j} = kron(speye(obj.dS),obj.bops{j});
            %             end
        end
        
        function createBopsFinLadder(obj)
            obj.bops{1} = sparse(diag(ones(1,obj.dB-1),-1));
            obj.bops{2} = sparse(obj.bops{1}');
            obj.bops{3} = sparse(diag(0:obj.dB-1));
            %             obj.Bops{3} = kron(speye(obj.dS),obj.bops{3});
        end
        
        function createBopsSpin(obj)
            J = (obj.dB-1)/2;
            m = (-J:J).';
            t = sqrt((J-m).*(J+m+1));
            obj.bops{1} = sparse(diag(t(1:end-1),-1));
            obj.bops{2} = sparse(obj.bops{1}');
            obj.bops{3} = sparse(diag(0:obj.dB-1)); %shift energy ground state = 0
        end
        
        function createOperators(obj)
            obj.createBops()
            obj.createSops();
            for i = 1:length(obj.bops)
                obj.Bops{i} = kron(obj.bops{i},speye(obj.dS^obj.nS));
            end
            obj.HS0 = sparse(zeros(obj.dB*obj.dS^obj.nS));
            for j = 1:obj.nS
                dBtoS = obj.dS^(j-1);
                dStoE = obj.dS^(obj.nS-j);
                for k = 1:length(obj.sops)
                    obj.Sops{k,j} = kron(speye(obj.dB*dBtoS),kron(obj.sops{k},speye(dStoE)));
                end
                obj.HS0 = obj.HS0 + obj.Sops{3,j};
            end
            obj.create2BodyExc();
        end
        
        function rho = createGaussianBatt(obj,mu,sigma)
            m = diag(obj.bops{3});
            dist = exp(-1/2/sigma^2*(m-mu).^2);
            dist = dist/sum(dist);
            rho = diag(dist);
        end
        
        function createInitState(obj,rhoS,rhoB)
            obj.rhoB{1} = sparse(rhoB);
            obj.rhoIndSth = sparse(rhoS);
            obj.rhoSth = rhoS;
            for j = 2:obj.nS
                obj.rhoSth = kron(obj.rhoSth,rhoS);
            end
            obj.rhoSth = sparse(obj.rhoSth);
            obj.rhoBS0 = sparse(kron(rhoB,obj.rhoSth));
            obj.EB = []; %reset energies
            obj.ErgB = [];
            obj.getBattEnergy();
        end
        
        function runSameInteraction(obj,Nstep,V,t)
            obj.U = expm(-1i*V*t);
            for j = 1:Nstep
                obj.runSingleStep();
            end
        end
        
        function runSameInteractionwReset(obj,Nstep,V,t)
            obj.U = expm(-1i*V*t);
            for j = 1:Nstep
                obj.runSingleStep();
                obj.rhoBS0 = kron(obj.rhoB{end},obj.rhoSth);
            end
        end
        
        function runDiffInteraction(obj,V,t)
            if (length(t) == 1)
                for j = 1:length(V)
                    obj.U = expm(-1i*V{j}*t);
                    obj.runSingleStep();
                end
            elseif (length(V)==1)
                for j = 1:length(t)
                    obj.U = expm(-1i*V*t(j));
                    obj.runSingleStep();
                end
            elseif (length(t)==length(V))
                for j = 1:length(t)
                    obj.U = expm(-1i*V{j}*t(j));
                    obj.runSingleStep();
                end
            else
                error('Dimensions of t and V don''t match');
            end
        end
        
        function runDiffInteractionwReset(obj,V,t)
            if (length(t) == 1)
                for j = 1:length(V)
                    obj.U = expm(-1i*V{j}*t);
                    obj.runSingleStep();
                    obj.rhoBS0 = kron(obj.rhoB{end},obj.rhoSth);
                end
            elseif (length(V)==1)
                for j = 1:length(t)
                    obj.U = expm(-1i*V*t(j));
                    obj.runSingleStep();
                    obj.rhoBS0 = kron(obj.rhoB{end},obj.rhoSth);
                end
            elseif (length(t)==length(V))
                for j = 1:length(t)
                    obj.U = expm(-1i*V{j}*t(j));
                    obj.runSingleStep();
                    obj.rhoBS0 = kron(obj.rhoB{end},obj.rhoSth);
                end
            else
                error('Dimensions of t and V don''t match');
            end
        end
        
        function create2BodyExc(obj)
            for j = 1:obj.nS
                obj.twoBodyV{j} = obj.Bops{1}*obj.Sops{2,j};
                obj.twoBodyV{j} = obj.twoBodyV{j} + obj.twoBodyV{j}';
            end
        end
        
        function createParallelExc(obj)
            obj.parV = sparse(zeros(obj.dB*obj.dS^obj.nS));
            for j = 1:obj.nS
                obj.parV = obj.parV + obj.twoBodyV{j};
            end
            obj.parV = obj.parV + obj.parV';
        end
        
        function runSingleStep(obj)
            obj.rhoBS0 = obj.U*obj.rhoBS0*obj.U';
            obj.getBattEnergy();
        end
        
        function getBattEnergy(obj)
            [obj.rhoB{end+1},~] = obj.ptrace(obj.rhoBS0,obj.dB,obj.dS^obj.nS);
            obj.EB = [obj.EB obj.getExpectation(obj.rhoB{end},obj.bops{3})];
            obj.ErgB = [obj.ErgB obj.getErgotropy(obj.rhoB{end},obj.bops{3})];
        end
        
        function getInvTEnergy(obj)
            obj.Eth = obj.getExpectation(obj.rhoIndSth,obj.sops{3});
            obj.Ergth = obj.getErgotropy(obj.rhoIndSth,obj.sops{3});
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
            %             [a,b] = size(H);
            %             [c,d] = size(rho);
            %             if (a~=b)||(c~=d)||(a~=c)
            %                 error('Check input matrices.');
            %             end
            E = real(full(sum(sum(rho.*(H.'))))); %more efficient than predefined tr(rho*H)
        end
        
        function resetState(obj)
            tmp = obj.rhoB{1};
            clear obj.rhoB
            obj.rhoB = cell(0);
            obj.rhoB{1} = tmp;
        end
        
        function E = getErgotropy(obj,rho,H)
            %getErgotropy computes the ergotropy of all states of rhoB with
            %respect to the battery bare Hamiltonian.
            rho = full(rho);
            [pvec, p] = eig(rho);
            [p,idx] = sort(diag(p),'descend');
            pvec = pvec(:,idx);
            [evec, e] = eig(full(H));
            E = real(obj.getExpectation(rho,H) - sum(p.*diag(e)));
            
        end
        
        
        
    end
    
end

