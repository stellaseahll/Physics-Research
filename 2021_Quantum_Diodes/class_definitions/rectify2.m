%SS 29.01.2021
%Class definition for rectifier consisting of N non-interacting qubits, one
%bath couples to them collectively, the other couples to them locally

classdef rectify2 < handle
    
    properties
        N = 2;
        dS = 4; %total dimension
        % all in units of qubit gap = 1
        k1 = 0.001; %bath coupling to spin 1
        k2 = 0.001; %bath coupling to spin 2
        T1 = 1; %bath 1 temperature
        T2 = 1; %bath 2 temperature
        %% General Operators
        pauliOps
        spinOps
        collOps
        %% Hamiltonians
        H0 %bare system Hamiltonian
        Hint %interaction Hamiltonian between 2 spins
        HS %total system Hamiltoninan
        
        %% Liouville Ops
        U
        dissLocP
        dissLocM
        dissCollP
        dissCollM
        
        %% energy operactors
        enLocP
        enLocM
        enCollP
        enCollM
        %% Bath operators
        Qop
        %% Total Liouvillian
        L %1: total Liouvillian of the master equation for T1-S1-S2-T2, 2:total Liouvillian of the master equation for T2-S1-S2-T1
        %% Steady state
        rhoSS %1: steady state for T1-S1...T2, 2: for  T2-S1...T1  - = COLLECTIVE ... = LOCAL
        Qexp 
    end
    
    methods
        function obj = rectify2(N,k1,k2,T1,T2)
            obj.N = N;
            obj.dS = 2^N;
            obj.k1 = k1;
            obj.k2 = k2;
            if (T2>T1)
                tmp = T2;
                T2 = T1;
                T1 = tmp;
                warning('Input temperature is swapped such that T1>T2.');
            end
            obj.T1 = T1;
            obj.T2 = T2;
            %%Create system operators
            obj.createSpinOps();
            obj.createOperators();
            %% Create Liouville
            obj.createLOps();
%             
            %% Find steady state
            obj.findSS();
            obj.getHeatCurrent();
        end
        
        function createSpinOps(obj)
            obj.pauliOps{1} = sparse([0,0;1,0]);
            obj.pauliOps{2} = sparse([0,1;0,0]);
            obj.pauliOps{3} = sparse([0 0;0 1]);
            obj.collOps{1} = sparse(zeros(2^obj.N));
            obj.collOps{2} = sparse(zeros(2^obj.N));
            obj.collOps{3} = sparse(zeros(2^obj.N));
            for j = 1:obj.N
                dStoj = 2^(j-1);
                djtoE = 2^(obj.N-j);
                for k = 1:3
                    obj.spinOps{k,j} = sparse(kron(speye(dStoj),kron(obj.pauliOps{k},speye(djtoE))));
                    obj.collOps{k} = obj.collOps{k} + obj.spinOps{k,j};
                end
            end
        end
        
        function createOperators(obj)
            %assume exchange interaction for now
            obj.HS = obj.collOps{3};
            obj.U = 1i*( obj.rightMultiply(obj.HS) - obj.leftMultiply(obj.HS));
            obj.dissLocP = obj.diss(obj.spinOps{1,1});
            obj.dissLocM = obj.diss(obj.spinOps{2,1});
            obj.enLocP = obj.dissEnergy(obj.spinOps{1,1},obj.HS);
            obj.enLocM = obj.dissEnergy(obj.spinOps{2,1},obj.HS);
            for j = 2:obj.N
                obj.dissLocP = obj.dissLocP + obj.diss(obj.spinOps{1,j});
                obj.dissLocM = obj.dissLocM + obj.diss(obj.spinOps{2,j});
                obj.enLocP = obj.enLocP +  obj.dissEnergy(obj.spinOps{1,j},obj.HS);
                obj.enLocM = obj.enLocM +  obj.dissEnergy(obj.spinOps{2,j},obj.HS);
            end
            obj.dissCollP = obj.diss(obj.collOps{1});
            obj.dissCollM = obj.diss(obj.collOps{2});
            obj.enCollP = obj.dissEnergy(obj.collOps{1},obj.HS);
            obj.enCollM = obj.dissEnergy(obj.collOps{2},obj.HS);
        end
        
        function createLOps(obj)

            %forward
            nf1 = 1/exp(1/obj.T1);
            nf2 = 1/exp(1/obj.T2);
            obj.L{1} = obj.U + obj.k1*nf1*obj.dissCollP + obj.k1*(nf1+1)*obj.dissCollM +...
                obj.k2*nf2*obj.dissLocP + obj.k2*(nf2+1)*obj.dissLocM;
            %backward
            nb1 = 1/exp(1/obj.T2);
            nb2 = 1/exp(1/obj.T1);
            obj.L{2} = obj.U + obj.k2*nb1*obj.dissCollP + obj.k2*(nb1+1)*obj.dissCollM +...
                obj.k1*nb2*obj.dissLocP + obj.k1*(nb2+1)*obj.dissLocM;
            %%get Qops
            obj.Qop{1} = obj.k1*nf1*obj.enCollP + obj.k1*(nf1+1)*obj.enCollM;
            obj.Qop{2} = obj.k1*nb2*obj.enLocP + obj.k1*(nb2+1)*obj.enLocM;
            
        end
        
        function findSS(obj) %finds the steady state for both Lf and Lb
            for j = 1:2
                SS = null(full(obj.L{j}));
                SS = reshape(SS,[obj.dS,obj.dS]);
                SS = SS/trace(SS);
                obj.rhoSS{j} = SS;
            end
        end
        
        function getHeatCurrent(obj)
            for j = 1:2
                obj.Qexp(j) = obj.getExpectation(obj.rhoSS{j},obj.Qop{j});
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
        
        function L = diss(obj,A)
            tmp = A'*A;
            L = obj.lrMultiply(A) - 0.5*obj.leftMultiply(tmp) - 0.5*obj.rightMultiply(tmp);
        end
        
        function Q = dissEnergy(obj,A,H)
            tmp = A'*A;
            Q = A'*H*A - 0.5*tmp*H - 0.5*H*tmp;
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
        
        
        
        
        
    end
    
end