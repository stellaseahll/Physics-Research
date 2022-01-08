% Class definition for N1 qubits interacting with N2 qubits.
% The N1 and N2 qubits are separately in thermal contact with baths of
% temperatures T1 and T2 with rates k1 and k2.

classdef spCollectiveSpins < handle
    
    properties
        %% system parameters
        N1 = 1; %default single qubit
        N2 = 1;
        w1 = 1;  %frequency
        w2 = 1;
        k1 = 0.001; %thermalization rate
        k2 = 0.001;
        T1 %temperature
        T2
        n1 %mean occupation
        n2
        dTot %total dimension of N1+N2
        %% spin matrices
        smat %pauli spin ops cell size 1x7
        J1mat %x,y,z,p,m,pm,mp collective spin matrix cell size: 1x7
        J2mat
        indJ1mat %x,y,z,p,m,pm,mp individual spin operator matrix cell size: N1x7
        indJ2mat %x,y,z,p,m,pm,mp
        %% Hamiltonian
        H0
        H
        indH
        U
        indU
        %% Dissipators
        L1 %collective bath system coupling
        Q1 %to get heat flow
        L2
        Q2
        indL1 %individual bath system coupling
        indQ1
        indL2
        indQ2
        %% Full Liouville
        Ltot
        indLtot
        %% Projectors (for multiple subspaces)
        P
        indP
        eigSS
        indeigSS
        %% States
        rho0
        rhoSS
        indrhoSS
    end
    
    methods
        function obj = spCollectiveSpins(N,T,k,w)
            obj.N1 = N(1);
            obj.N2 = N(2);
            obj.dTot = 2^(sum(N));
            obj.T1 = T(1);
            obj.T2 = T(2);
            if (nargin>=3)
                obj.k1 = k(1);
                obj.k2 = k(2);
            end
            if (nargin>=4)
                obj.w1 = w(1);
                obj.w2 = w(2);
            end
            obj.n1 = 1/(exp(obj.w1/obj.T1)-1);
            obj.n2 = 1/(exp(obj.w2/obj.T2)-1);
            obj.getPauliOperators();
            obj.getIndividualOperators();
            obj.getCollectiveOperators();
            obj.getBareHamiltonian();
            obj.getCollectiveDiss();
            obj.getIndividualDiss();
            
        end
        
        function getPauliOperators(obj)
            obj.smat{1} = sparse([0 1; 1 0]);
            obj.smat{2} = sparse([0 1i; -1i 0]);
            obj.smat{3} = sparse([-0.5 0 ; 0 0.5]);
            obj.smat{4} = sparse([0 0; 1 0]);
            obj.smat{5} = sparse([0 1; 0 0]);
            obj.smat{6} = sparse([0 0; 0 1]);
            obj.smat{7} = sparse([1 0; 0 0]);
        end
        
        function getIndividualOperators(obj)
            %system 1
            for j = 1:obj.N1
                dToptoCurr = 2^(j-1);
                dCurrtoEnd = 2^(obj.N1+obj.N2-j);
                for k = 1:length(obj.smat)
                    obj.indJ1mat{j,k} = kron(speye(dToptoCurr),kron(obj.smat{k},speye(dCurrtoEnd)));
                end
            end
            %system 2
            for j = 1:obj.N2
                dToptoCurr = 2^(obj.N1+j-1);
                dCurrtoEnd = 2^(obj.N2-j);
                for k = 1:length(obj.smat)
                    obj.indJ2mat{j,k} = kron(speye(dToptoCurr),kron(obj.smat{k},speye(dCurrtoEnd)));
                end
            end
        end
        
        function getCollectiveOperators(obj)
            %system 1
            for k = 1:5
                obj.J1mat{k} = obj.indJ1mat{1,k};
                for j = 2:obj.N1
                    obj.J1mat{k} = obj.J1mat{k} + obj.indJ1mat{j,k};
                end
            end
            obj.J1mat{6} = obj.J1mat{4}*obj.J1mat{5};
            obj.J1mat{7} = obj.J1mat{5}*obj.J1mat{4};
            %system 2
            for k = 1:5
                obj.J2mat{k} = obj.indJ2mat{1,k};
                for j = 2:obj.N2
                    obj.J2mat{k} = obj.J2mat{k} + obj.indJ2mat{j,k};
                end
            end
            obj.J2mat{6} = obj.J2mat{4}*obj.J2mat{5}; %J+*J- \neq spm1+spm2...
            obj.J2mat{7} = obj.J2mat{5}*obj.J2mat{4};
        end
        
        function getBareHamiltonian(obj)
            obj.H0 = obj.J1mat{3}*obj.w1 + obj.J2mat{3}*obj.w2;
            obj.H = obj.H0;
            obj.indH = obj.H0;
        end
        
        function getCollectiveDiss(obj)
            obj.L1 = obj.k1*(obj.n1+1)*(obj.lrMultiply(obj.J1mat{5}) - ...
                0.5*( obj.leftMultiply(obj.J1mat{6}) + obj.rightMultiply(obj.J1mat{6})))+...
                obj.k1*obj.n1*(obj.lrMultiply(obj.J1mat{4}) - ...
                0.5*( obj.leftMultiply(obj.J1mat{7}) + obj.rightMultiply(obj.J1mat{7})));
            obj.L2 = obj.k2*(obj.n2+1)*(obj.lrMultiply(obj.J2mat{5}) - ...
                0.5*( obj.leftMultiply(obj.J2mat{6}) + obj.rightMultiply(obj.J2mat{6})))+...
                obj.k2*obj.n2*(obj.lrMultiply(obj.J2mat{4}) - ...
                0.5*( obj.leftMultiply(obj.J2mat{7}) + obj.rightMultiply(obj.J2mat{7})));
        end
        
        function getIndividualDiss(obj)
            %system 1
            obj.indL1 = zeros(obj.dTot^2);
            for j = 1:obj.N1
                obj.indL1 = obj.indL1 + obj.k1*(obj.n1+1)*(obj.lrMultiply(obj.indJ1mat{j,5}) - ...
                    0.5*( obj.leftMultiply(obj.indJ1mat{j,6}) + obj.rightMultiply(obj.indJ1mat{j,6})))+...
                    obj.k1*obj.n1*(obj.lrMultiply(obj.indJ1mat{j,4}) - ...
                    0.5*( obj.leftMultiply(obj.indJ1mat{j,7}) + obj.rightMultiply(obj.indJ1mat{j,7})));
            end
            %system 2
            obj.indL2 = zeros(obj.dTot^2);
            for j = 1:obj.N2
                obj.indL2 = obj.indL2 + obj.k2*(obj.n2+1)*(obj.lrMultiply(obj.indJ2mat{j,5}) - ...
                    0.5*( obj.leftMultiply(obj.indJ2mat{j,6}) + obj.rightMultiply(obj.indJ2mat{j,6})))+...
                    obj.k2*obj.n2*(obj.lrMultiply(obj.indJ2mat{j,4}) - ...
                    0.5*( obj.leftMultiply(obj.indJ2mat{j,7}) + obj.rightMultiply(obj.indJ2mat{j,7})));
            end
        end
        
        function getCollectiveU(obj)
            obj.U = 1i*(obj.rightMultiply(obj.H) - obj.leftMultiply(obj.H));
        end
        
        function getIndividualU(obj)
            obj.indU = 1i*(obj.rightMultiply(obj.indH) - obj.leftMultiply(obj.indH));
        end
        
        function getCollectiveLiouville(obj)
            obj.getCollectiveU();
            obj.Ltot = obj.U + obj.L1 + obj.L2;
        end
        
        function getIndividualLiouville(obj)
            obj.getIndividualU();
            obj.indLtot = obj.indU + obj.indL1 + obj.indL2;
        end
        
        function M = getQOp(obj,k,op)
            tmp = op'*op ;
            M = k*(op'*obj.H*op - 0.5*obj.H*tmp - 0.5*tmp*obj.H);
        end
        
        function getCollectiveQOps(obj)
            obj.Q1 = obj.getQOp(obj.k1*(obj.n1+1),obj.J1mat{5})+obj.getQOp(obj.k1*obj.n1,obj.J1mat{4});
            obj.Q2 = obj.getQOp(obj.k2*(obj.n2+1),obj.J2mat{5})+obj.getQOp(obj.k2*obj.n2,obj.J2mat{4});    
        end
        
        function getIndividualQOps(obj)
            obj.indQ1 = zeros(obj.dTot);
            obj.indQ2 = zeros(obj.dTot);
            for j = 1:obj.N2
                obj.indQ1 = obj.indQ1 + obj.getQOp(obj.k1*(obj.n1+1),obj.indJ1mat{j,5})+obj.getQOp(obj.k1*obj.n1,obj.indJ1mat{j,4});
                obj.indQ2 = obj.indQ2 + obj.getQOp(obj.k2*(obj.n2+1),obj.indJ2mat{j,5})+obj.getQOp(obj.k2*obj.n2,obj.indJ2mat{j,4});
            end
        end
        
        function getAllQOps(obj)
            obj.getIndividualQOps();
            obj.getCollectiveQOps();
        end
        
        function [obst p0] = runME(obj,dt,obs)
            % Solves for state at specified time t using Liouvillian L
            % given initial state
            % Need to get L before running this
            Nt = length(t);
            Nobs = length(obs);
            obst{1} = zeros(Nobs,Nt);
            for n=1:Nobs
                obst{1}(n,1) = obj.getExpectation(obj.rho0,obs{n});
            end
            obst{2} = obst{1};
            rho0 = reshape(obj.rho0,obj.dTot^2,1);
            
            for j = 2:Nt
                Rhot = reshape(expm(obj.Ltot*t(j))*rho0,obj.dTot,obj.dTot);
                indRhot = reshape(expm(obj.indLtot*t(j))*rho0,obj.dTot,obj.dTot);
                for n=1:Nobs
                    obst{1}(n,j) = obj.getExpectation(Rhot,obs{n});
                    obst{2}(n,j) = obj.getExpectation(indRhot,obs{n});
                end
            end
        end
        
        function getCollectiveSSProjector(obj)
            eigS = spnull(obj.Ltot);
            [r,c] = size(eigS);
            for j = 1:c
                rho = reshape(eigS(:,j),obj.dTot,obj.dTot);
                absSS(j) = max(max(abs(rho-rho')));
            end
            obj.P = zeros(size(obj.Ltot));
            idx = find(absSS<1e-12);
            for j = idx
                obj.P = obj.P + eigS(:,j)*eigS(:,j)';
            end
            obj.eigSS = eigS(:,idx);
            
        end
        
        function getIndividualSSProjector(obj)
            eigS = spnull(obj.indLtot);
            [r,c] = size(eigS);
            for j = 1:c
                rho = reshape(eigS(:,j),obj.dTot,obj.dTot);
                absSS(j) = max(max(abs(rho-rho')));
            end
            obj.indP = zeros(size(obj.indLtot));
            idx = find(absSS<1e-12);
            for j = idx
                obj.indP = obj.indP + eigS(:,j)*eigS(:,j)';
            end
            obj.indeigSS = eigS(:,idx);
        end
        
        function getAllProjectors(obj)
            obj.getCollectiveSSProjector();
            obj.getIndividualSSProjector();
        end
        
        
        function getCollectiveSS(obj)
            tmp = obj.P*reshape(obj.rho0,obj.dTot^2,1);
            tmp = reshape(tmp,obj.dTot,obj.dTot);
            tmp = tmp/trace(tmp);
            obj.rhoSS = (tmp+tmp')/2;
            
        end
        
        function getIndividualSS(obj)
            tmp = obj.indP*reshape(obj.rho0,obj.dTot^2,1);
            tmp = reshape(tmp,obj.dTot,obj.dTot);
            tmp = tmp/trace(tmp);
            obj.indrhoSS = (tmp+tmp')/2;
        end
        
        function getAllSS(obj)
             obj.getIndividualSS();
             obj.getCollectiveSS();
        end
        
        function obst=computeSSvalues(obj,obs)
            obj.getAllQOps();
            obj.getAllSS();
            Nobs = length(obs);
            obst = zeros(Nobs+2,2);
            for n=1:Nobs
                obst(n,1) = obj.getExpectation(obj.rhoSS,obs{n});
                obst(n,2) = obj.getExpectation(obj.indrhoSS,obs{n});
            end
            obst(Nobs+1,1) = obj.getExpectation(obj.rhoSS,obj.Q1);
            obst(Nobs+2,1) = obj.getExpectation(obj.rhoSS,obj.Q2);
            obst(Nobs+1,2) = obj.getExpectation(obj.indrhoSS,obj.indQ1);
            obst(Nobs+2,2) = obj.getExpectation(obj.indrhoSS,obj.indQ2);
        end
        
        function M = leftMultiply(obj,A)
            M = (kron(eye(obj.dTot),A));
        end
        
        function M = rightMultiply(obj,A)
            M = (kron(A.',eye(obj.dTot)));
        end
        
        function M = lrMultiply(obj,A)
            M = (kron(conj(A),A));
        end
        
        function E = getExpectation(obj,rho,H)
            %getExpectation computes the expectation value tr(rho*H)
            [a,b] = size(H);
            [c,d] = size(rho);
            if (a~=b)||(c~=d)||(a~=c)
                error('Check input matrices.');
            end
            E = full(sum(sum(rho.*(H.')))); %more efficient than predefined tr(rho*H)
        end
    end
    
end


