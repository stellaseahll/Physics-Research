%SS 29.01.2021
%Class definition for rectifier consisting of 2 spins as system.

classdef rectify < handle
    
    properties
        dS1 = 2; %spin 1 dimension
        dS2 = 2; %spin 2 dimension
        dS = 4; %total dimension
        omega1 = 1; %spin 1 frequency
        omega2 = 1; %spin 2 frequency
        k1 = 0.001; %bath coupling to spin 1
        k2 = 0.001; %bath coupling to spin 2
        T1 = 1; %bath 1 temperature
        T2 = 1; %bath 2 temperature
        g = 0.01; %spin spin coupling between 2 spins
        %% General Operators
        S1ops
        S2ops
        
        %% Hamiltonians
        H0 %bare system Hamiltonian
        Hint %interaction Hamiltonian between 2 spins
        HS %total system Hamiltoninan
        
        %% Bath operators
        Qop
        %% Total Liouvillian
        L %1: total Liouvillian of the master equation for T1-S1-S2-T2, 2:total Liouvillian of the master equation for T2-S1-S2-T1
        %% Steady state
        rhoSS %1: steady state for T1-S1-S2-T2, 2: for  T2-S1-S2-T1
        Qexp
    end
    
    methods
        function obj = rectify(ds1,ds2,om1,om2,g,k1,k2,T1,T2)
            obj.dS1 = ds1;
            obj.dS2 = ds2;
            obj.dS = ds1*ds2;
            obj.omega1 = om1;
            obj.omega2 = om2;
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
            obj.g = g;
            %%Create system operators
            obj.createSpinOps1();
            obj.createSpinOps2();
            
            %% Create Liouville
            obj.createHamiltonian();
            obj.createLOps();
            
            %% Find steady state
            obj.findSS();
            obj.getHeatCurrent();
        end
        
        function createSpinOps1(obj)
            % create operators of S1 and S2
            % index: 1 = p, 2 = m, 3 = z
            J = (obj.dS1-1)/2;
            m = (-J:J).';
            t = sqrt((J-m).*(J+m+1));
            obj.S1ops{1} = sparse(diag(t(1:end-1),-1));
            obj.S1ops{2} = obj.S1ops{1}';
            obj.S1ops{3} = sparse(diag(m));
            for j = 1:3
                obj.S1ops{j} = kron(obj.S1ops{j},speye(obj.dS2));
            end
        end
        
        function createSpinOps2(obj)
            % create operators of S1 and S2
            % index: 1 = p, 2 = m, 3 = z
            J = (obj.dS2-1)/2;
            m = (-J:J).';
            t = sqrt((J-m).*(J+m+1));
            obj.S2ops{1} = sparse(diag(t(1:end-1),-1));
            obj.S2ops{2} = obj.S2ops{1}';
            obj.S2ops{3} = sparse(diag(m));
            for j = 1:3
                obj.S2ops{j} = kron(speye(obj.dS1),obj.S2ops{j});
            end
        end
        
        function createHamiltonian(obj)
            %assume exchange interaction for now
            obj.H0 = obj.omega1*obj.S1ops{3} + obj.omega2*obj.S2ops{3};
            obj.Hint = obj.g*obj.S1ops{1}*obj.S2ops{2};
            obj.Hint = obj.Hint + obj.Hint';
            obj.HS = obj.H0 + obj.Hint;
        end
        
        function createLOps(obj)
            U = 1i*( obj.rightMultiply(obj.HS) - obj.leftMultiply(obj.HS));
            %forward
            nf1 = 1/exp(obj.omega1/obj.T1);
            nf2 = 1/exp(obj.omega2/obj.T2);
            obj.L{1} = U + obj.k1*nf1*obj.diss(obj.S1ops{1}) + obj.k1*(nf1+1)*obj.diss(obj.S1ops{2}) +...
                obj.k2*nf2*obj.diss(obj.S2ops{1}) + obj.k2*(nf2+1)*obj.diss(obj.S2ops{2});
            %backward
            nb1 = 1/exp(obj.omega1/obj.T2);
            nb2 = 1/exp(obj.omega2/obj.T1);
            obj.L{2} = U + obj.k2*nb1*obj.diss(obj.S1ops{1}) + obj.k2*(nb1+1)*obj.diss(obj.S1ops{2}) +...
                obj.k1*nb2*obj.diss(obj.S2ops{1}) + obj.k1*(nb2+1)*obj.diss(obj.S2ops{2});
            %%get Qops
            obj.Qop{1} = obj.k1*nf1*obj.dissEnergy(obj.S1ops{1},obj.HS) + obj.k1*(nf1+1)*obj.dissEnergy(obj.S1ops{2},obj.HS);
            obj.Qop{2} = obj.k1*nb2*obj.dissEnergy(obj.S2ops{1},obj.HS) + obj.k1*(nb2+1)*obj.dissEnergy(obj.S2ops{2},obj.HS);
            
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