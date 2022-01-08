
classdef model3SpinsFridgeTwoPhotonHW < handle
    % gamma matrix separate from expectation values
    % For all combinations
    % Hot, Cold, Work
    properties
        %eigenfrequencies of qubits
        E1=5;
        E2=1;
        E3=4;
        g=0.1; %3-spin coupling rate (only resonant RWA terms!)
        % temperature terms kB*Ti
        T1 = 2.0;
        T2 = 1.0; 
        T3 = 5.0;
        
        %bath coupling rates
        k1 = 1e-3;
        k2 = 1e-3;
        k3 = 1e-3;
        
        %spin number
        J = 0.5;
        dim = 2;
        %operators
        LOps = cell(0); %cell array of Lindblad operators
        coeffLOps = [];
        bathIdx; %to keep track of different bath operators
        H %Hamilton operator
        Hint % interaction Hamilton operator
        %H0 % interaction Hamilton operator
        U %Eigen vectors of H
        E % eigenvalues of H
        L %total Liouvillian of the master equation
        P
        eigVal = []; %eigenvalues of gamma matrix
        
        %bath type
        bathType = 'flat'; %bath type: 'flat' or 'ohmic'
        intType = 'RWA';  %interaction type: 'RWA' or 'full'
        meType = 'partial'; %master equation type: 'partial','local','global','breuer'
        %states
        rho % steady state, i.e. kernel element of L
        rhoH
        rhoC
        rhoW
        rhoHC
        rhoCW
        rhoHW
        rhoSS
        %permutation
        %coarsegrain time
        deltat 
        gamma = cell(1,3)
    end
    
    methods
        
        function obj = model3SpinsFridgeTwoPhotonHW(G,Ts,ks,Es,Btype,Itype,Mtype,J)
            if nargin>=1
                obj.g = G;
            end
            if nargin>=2
                obj.T1 = Ts(1);
                obj.T2 = Ts(2);
                obj.T3 = Ts(3);
            end
            if nargin>=3
                obj.k1 = ks(1);
                obj.k2 = ks(2);
                obj.k3 = ks(3);
            end
            if nargin>=4
                obj.E1 = Es(1);
                obj.E2 = Es(2);
                obj.E3 = Es(3);
            end
            if nargin>=5
                switch Btype
                    case {'flat','f','F'}
                        obj.bathType = 'flat';
                    case {'ohmic','ohm','o','O'}
                        obj.bathType = 'ohmic';
                    case {'cubic','cube','c','C'}
                        obj.bathType = 'cubic';
                end
            end
            if nargin>=6
                switch Itype
                    case {'RWA','rwa','r','R'}
                        obj.intType = 'RWA';
                    case {'full','f','F','x','X','xxx','XXX'}
                        obj.intType = 'full';
                end
            end
            if nargin>=7
                switch Mtype
                    case {'partial','P','p'}
                        obj.meType = 'partial';
                    case {'local','L','l'}
                        obj.meType = 'local';
                    case {'global','G','g'}
                        obj.meType = 'global';
                    case {'breuer','B','b'}
                        obj.meType = 'breuer';
                end
            end
            if nargin>=8
                obj.J = J;
                obj.dim = 2*J+1;
            end
            obj.deltat = 1;%./sqrt(obj.g*min([obj.E1,obj.E2,obj.E3]));
            
            obj.resetOperators();
        end
%       
%         function tmp = getOp(obj,eigvec,i,j,X)
%             tmp = (eigvec(:,i)*eigvec(:,j)') * (eigvec(:,i)'*X*eigvec(:,j));
%         end

        % check definition for sigmax
        function sp = getSigmap(obj)
            m = -obj.J:obj.J;
            t = sqrt((obj.J-m).*(obj.J+m+1));
            sp = diag(t(1:end-1),-1)/sqrt(2);
        end
        
        function sm = getSigmam(obj)
            m = -obj.J:obj.J;
            t = sqrt((obj.J+m).*(obj.J-m+1));
            sm = diag(t(2:end),1)/sqrt(2);
        end
        
        
        function resetOperators(obj) %create all operators and superops.
            s_p = obj.getSigmap();
            s_m = obj.getSigmam();
            s_z = diag(-obj.J:obj.J);
            s_x = (s_p+s_m);
            Id = eye(2*obj.J+1);

%             n1 = nbar(obj,obj.E1,obj.E2);
%             n2 = nbar(obj,obj.E2,obj.T2);
%             n3 = nbar(obj,obj.E3,obj.T3); 

            switch obj.intType
                case 'RWA'
                    obj.Hint = obj.g*kron( s_p, kron(s_m, s_m) );
                    obj.Hint = obj.Hint + obj.Hint';
                case 'full'
                    obj.Hint = obj.g*kron( s_x, kron(s_x, s_x) );
            end
            obj.H = obj.Hint + obj.E1*kron( s_z, kron(Id, Id) ) + ...
                    obj.E2*kron( Id, kron(s_z, Id) ) + ...
                    obj.E3*kron( Id, kron(Id, s_z) );
            [obj.U , obj.E] = eig(obj.H,'vector');
            obj.U(abs(obj.U)<10*eps) = 0;
            
            obj.U = sparse(obj.U);
            obj.L = -1i*leftMultiply(obj.H,obj.dim^3) + 1i*rightMultiply(obj.H,obj.dim^3);
            %obj.LH = -1i*leftMultiply(obj.H,8) + 1i*rightMultiply(obj.H,8);
            switch obj.meType
                case 'local'
                    obj.bathIdx = [0 2 4 6];
                    n1 = nbar(obj,obj.E1*2,obj.T1);
                    n2 = nbar(obj,obj.E2,obj.T2);
                    n3 = nbar(obj,obj.E3*2,obj.T3);
                    obj.LOps{1} = sqrt(obj.k1*2*(n1+1)) * kron( s_m*s_m, kron(Id, Id) );
                    obj.LOps{2} = sqrt(obj.k1*2*n1) * kron( s_p*s_p, kron(Id, Id) );
                    obj.LOps{3} = sqrt(obj.k2*(n2+1)) * kron( Id, kron(s_m, Id) );
                    obj.LOps{4} = sqrt(obj.k2*n2) * kron( Id, kron(s_p, Id) );
                    obj.LOps{5} = sqrt(obj.k3*2*(n3+1)) * kron( Id, kron(Id, s_m*s_m) );
                    obj.LOps{6} = sqrt(obj.k3*2*n3) * kron( Id, kron(Id, s_p*s_p) );
                otherwise
                    obj.bathIdx = 0;
                    switch obj.bathType
                        case 'cubic'
                            getLindblad(obj,kron(s_x,kron(Id,Id)),obj.k1/(obj.E1)^3,obj.T1);
                            getLindblad(obj,kron(Id,kron(s_x,Id)),obj.k2/(obj.E2)^3,obj.T2);
                            getLindblad(obj,kron(Id,kron(Id,s_x)),obj.k3/(obj.E3)^3,obj.T3);
                        case 'ohmic'
                            getLindblad(obj,kron(s_x,kron(Id,Id)),obj.k1/obj.E1,obj.T1);
                            getLindblad(obj,kron(Id,kron(s_x,Id)),obj.k2/obj.E2,obj.T2);
                            getLindblad(obj,kron(Id,kron(Id,s_x)),obj.k3/obj.E3,obj.T3);
                        case 'flat'
                            getLindblad(obj,kron(s_x,kron(Id,Id)),obj.k1,obj.T1);
                            getLindblad(obj,kron(Id,kron(s_x,Id)),obj.k2,obj.T2);
                            getLindblad(obj,kron(Id,kron(Id,s_x)),obj.k3,obj.T3);
                    end
            end
            %Hamiltonian part
            
            %add the Lindblad terms
            for j = 2:length(obj.bathIdx)   
                a=obj.bathIdx(j-1)+1;
                b=obj.bathIdx(j);
                for k = a:b
                    LL = sparse(obj.LOps{k});
                    tmp = LL'*LL;
                    obj.L = obj.L + rightMultiply(LL',obj.dim^3)*leftMultiply(LL,obj.dim^3) - ...
                    0.5*( leftMultiply(tmp,obj.dim^3) + rightMultiply(tmp,obj.dim^3));
                end
            end
            tmp = null(full(obj.L));
            [r,c] = size(tmp);
            obj.P = zeros(size(obj.L));
            for j = 1:c
                obj.P = obj.P + tmp(:,j)*tmp(:,j)';
            end
        end
        
        function Q = evalQCvaryGChi(obj,G,chi) 
            %evals -QC at varying g and chi (=obj.k2)
            %Does not save the SS als obj.rho!
            obj.g = G;
            chi0 = obj.k2;
            obj.k2 = chi;
            obj.k1 = chi/chi0*obj.k1;
            obj.k3 = chi/chi0*obj.k3;
            obj.resetOperators();
            SS = null(obj.L);
            SS = reshape(SS,[obj.dim^3,obj.dim^3]);
            SS = SS/trace(SS);
            a = obj.bathIdx(2)+1;
            b = obj.bathIdx(3);
            Q = 0;
            for k = a:b
                LL = obj.LOps{k};
                Q = Q + trace( obj.H * ( LL*(SS*(LL')) - 0.5*(LL'*LL)*SS - 0.5*SS*(LL'*LL) ) );
            end
            Q = -real(Q);
        end
        
        function [gMax,chiMax,QCmax,exitflag,output] = findMaxQCvaryGChi(obj)
            %constrained max searcher of QC as function of g and chi
            % chi is assumed to be equal to kappaC (i.e. omegaC=1)
            % initial guess is current g, chi set in obj.
            % constraints are: 1e-4 <= g <= 1, 1e-5 <= chi <= 0.1
            [xM,QCmax,exitflag,output] = fmincon(@(x)evalQCvaryGChi(obj,x(1),x(2)),[obj.g obj.k2], ...
                [],[],[],[],[1e-4, 1e-5],[1, 0.1]);
            gMax=xM(1);
            chiMax=xM(2);
            QCmax = -QCmax;
        end

        function findSS(obj,rhoIn) %finds the steady state of L
                tmp = obj.P*reshape(rhoIn,numel(rhoIn),1);
                tmp = reshape(tmp,length(rhoIn),length(rhoIn));
                obj.rho = tmp/trace(tmp);
%             [~, SS] = spspaces(sparse(obj.L),2,1e-10);
%             Q = SS{1};
%             J = SS{3};
%             SS = Q(:,J);
%             SS = reshape(full(SS),[obj.dim^3,obj.dim^3]);
%             SS = SS/trace(SS);
%             obj.rho = SS;
%             permMat = getPermMat(obj);
%             [obj.rhoH, obj.rhoCW] = ptrace(SS,obj.dim,obj.dim^2);
%             [obj.rhoHC, obj.rhoW] = ptrace(SS,obj.dim^2,obj.dim);
%             [obj.rhoHW, obj.rhoC] = ptrace(permMat*SS*permMat,obj.dim^2,obj.dim);
        end
        
        
        function permMat = getPermMat(obj)
            permMat = eye(obj.dim^3);
            for a = 1:obj.dim
                for b = 1:obj.dim
                     for c = b+1:obj.dim
                         id1 = ((a-1)*obj.dim+b-1)*obj.dim + c;
                         id2 = ((a-1)*obj.dim+c-1)*obj.dim + b;
                         permMat(id1,id1) = 0;
                         permMat(id2,id2) = 0;
                         permMat(id1,id2) = 1;
                         permMat(id2,id1) = 1;
                     end
                     
                end
            end
            permMat = sparse(permMat);
            
        end
        function rho = generateRhoIn(obj,E,T) %generate product of thermal states
            rho1 = exp(-E(1)/T(1)*(0:obj.dim-1));
            rho1 = diag(rho1/sum(rho1));
            rho2 = exp(-E(2)/T(2)*(0:obj.dim-1));
            rho2 = diag(rho2/sum(rho2));
            rho3 = exp(-E(3)/T(3)*(0:obj.dim-1));
            rho3 = diag(rho3/sum(rho3));
            rho = kron(rho1,kron(rho2,rho3));
        end
        
%         function rho = getThermalState(obj,bath) %reduced thermal state of spin 'bath'
%             switch bath
%                 case 1
%                     rho = diag([1 exp(-obj.E1/obj.T1)])/sum([1 exp(-obj.E1/obj.T1)]);
%                 case 2
%                     rho = diag([1 exp(-obj.E2/obj.T2)])/sum([1 exp(-obj.E2/obj.T2)]);
%                 case 3
%                     rho = diag([1 exp(-obj.E3/obj.T3)])/sum([1 exp(-obj.E3/obj.T3)]);
%             end
%         end
        
        
        function U=generateU(obj,dt)
            U = expm(obj.L*dt);
        end
        
        function y = timeEvolve(obj,rhoIn,dt,tmax)
            timestep = 0:dt:tmax;
            U=generateU(obj,dt);
            obj.rho = rhoIn;
            rhoVec = reshape(rhoIn,[64,1]);
            getReducedState1(obj);
            getReducedState2(obj);
            getReducedState3(obj);
            eigVal = eig(obj.rho,'vector');
            data = [getTotalHeatFlow(obj,1);getTotalHeatFlow(obj,2);getTotalHeatFlow(obj,3);...
                getEnergy(obj,1);getEnergy(obj,2);getEnergy(obj,3);sum(eigVal);min(eigVal)];
            for i = 1:length(timestep)-1
                rhoVec = U*rhoVec;
                obj.rho = reshape(rhoVec,[obj.dim^3 obj.dim^3]);
                getReducedState1(obj);
                getReducedState2(obj);
                getReducedState3(obj);
                eigVal = eig(obj.rho,'vector');
                tmp = [getTotalHeatFlow(obj,1);getTotalHeatFlow(obj,2);getTotalHeatFlow(obj,3);...
                    getEnergy(obj,1);getEnergy(obj,2);getEnergy(obj,3);sum(eigVal);min(eigVal)];
                data = [data tmp];
            end
            y = [timestep; data];
        end
        
        function n = nbar(obj,omega,T)
            n = 1./(exp(omega./T)-1);
            %cutoff criterion
%             n(T>100*omega) = T./omega(T>100*omega);
        end
        
        function Q = getHeatFlow(obj,j) %average energy change at SS coming from Lindblad term j
            %CHECK, modify as needed
            pi = obj.rho;
            LL = obj.LOps{j};
            Q = trace( obj.H * ( LL*(pi*(LL')) - 0.5*(LL'*LL)*pi - 0.5*pi*(LL'*LL) ) );
            
        end
        
        function Q = getTotalHeatFlow(obj,bath)
            
            a = obj.bathIdx(bath)+1;
            b = obj.bathIdx(bath+1);
            Q = 0;
            for k = a:b
                Q = Q + getHeatFlow(obj,k);
            end
            Q = real(Q);
        end
        
        function S = getEntropy(obj,rho)
            ev = real(eig(rho));
            S = -sum(ev.*log(ev));
        end
        
        function C = getCoherence(obj,rho)
            C = getEntropy(obj,diag(diag(rho)))- getEntropy(obj,rho);
        end
        
        function Cdata = getAllCoherence(obj)
            Ct = getCoherence(obj,obj.rho);
            Cl = getCoherence(obj,kron(obj.rhoH,kron(obj.rhoC,obj.rhoW)));
            Cg = Ct - Cl;
            Chc = getCoherence(obj,obj.rhoHC)-getCoherence(obj,kron(obj.rhoH,obj.rhoC));
            Chw = getCoherence(obj,obj.rhoHW)-getCoherence(obj,kron(obj.rhoH,obj.rhoW));
            Ccw = getCoherence(obj,obj.rhoCW)-getCoherence(obj,kron(obj.rhoC,obj.rhoW));
            Cbg = Chc+Chw+Ccw;
            Ch_cw = Ct-getCoherence(obj,kron(obj.rhoH,obj.rhoCW));
            Ctg = Ccw+Ch_cw;
            rho2 = obj.U'*obj.rho*obj.U;
            Cte = getCoherence(obj,rho2);
            Ct2 = getCoherence(obj,obj.rho(4:5,4:5));
            Ct2e = getCoherence(obj,rho2(4:5,4:5));
            Cdata = [Ct2 Ct2e Cte Ct Cg Cl Ctg Cbg];
            
        end
        
        function E = getEntanglement(obj)
            %WRONG! For correlated mixed states, E can be negative without
            %entanglement!
            warning('Entropy difference does not measure entanglement for mixed states!')
            E = getEntropy(obj,obj.rho)-getEntropy(obj,obj.rhoH)-getEntropy(obj,obj.rhoCW);
        end
        
        function [neg,logneg] = getNegativityHot(obj)
            %Measures entanglement by PT-negativity. neg > 0 means there
            %are negative eigenvalues of the PT of rho, i.e. entanglement
            %Use bipartition H - CW (specified my resonant interaction)
            [~,neg,logneg] = partialTransposeHot(obj.rho);
        end
        
        function [neg,logneg] = getNegativityCold(obj)
            %Measures entanglement by PT-negativity. neg > 0 means there
            %are negative eigenvalues of the PT of rho, i.e. entanglement
            %Use bipartition C - HW (specified by virtual qubit model)
            [~,neg,logneg] = partialTransposeCold(obj.rho);
        end
        
        function C = getConcurrences(obj)
            %Measures entanglement by concurrence > 0 for all biseparations
            %(following Brunner et al PRE 89, 2014)
            %entanglement only if positive.
            %TODO
        end

        function W = getGHZWitness(obj)
            %Measures entanglement by computing several GHZ witnesses,
            %negative values would indicate the GHZ-type entanglement.
            W.f000p111 = 3/4 - (obj.rho(1,1) + obj.rho(8,8) + obj.rho(1,8) + obj.rho(8,1))/2;
            W.f000pi111 = 3/4 - (obj.rho(1,1) + obj.rho(8,8) + 1i*obj.rho(1,8) - 1i*obj.rho(8,1))/2;
            W.f000m111 = 3/4 - (obj.rho(1,1) + obj.rho(8,8) - obj.rho(1,8) - obj.rho(8,1))/2;
            W.f100p011 = 3/4 - (obj.rho(5,5) + obj.rho(4,4) + obj.rho(5,4) + obj.rho(4,5))/2;
            W.f100pi011 = 3/4 - (obj.rho(5,5) + obj.rho(4,4) + 1i*obj.rho(5,4) - 1i*obj.rho(4,5))/2;
            W.f100m011 = 3/4 - (obj.rho(5,5) + obj.rho(4,4) - obj.rho(5,4) - obj.rho(4,5))/2;
            W.f110p001 = 3/4 - (obj.rho(7,7) + obj.rho(2,2) + obj.rho(7,2) + obj.rho(2,7))/2;
            W.f110pi001 = 3/4 - (obj.rho(7,7) + obj.rho(2,2) + 1i*obj.rho(7,2) - 1i*obj.rho(2,7))/2;
            W.f110m001 = 3/4 - (obj.rho(7,7) + obj.rho(2,2) - obj.rho(7,2) - obj.rho(2,7))/2;
        end
        
        
        function [A, B] = ptrace(AB, dimA, dimB)

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

        
        function E = getEnergy(obj,bath) %avg single-mode energy of spin 'bath' in SS 
            switch bath
                case 1
                    E = obj.rhoH(2,2)*obj.E1;
                case 2
                    E = obj.rhoC(2,2)*obj.E2;
                case 3
                    E = obj.rhoW(2,2)*obj.E3;
            end
            E = real(E);
            
        end
        
        function [trRho, minEig, maxImag] = testState(obj,numTrial)
            dt = logspace(-3,3,10);
            for i = 1:length(dt)
                U{i} = generateU(obj,dt(i));
            end
            trRho = zeros(numTrial,10);
            minEig = zeros(numTrial,10);
            maxImag = zeros(numTrial,10);

            for i = 1:numTrial
                randNum = rand(1,12)-0.5;
                psi = [randNum(1) + randNum(2)*1i; randNum(3)+randNum(4)*1i];
                psi = psi/norm(psi);
                A = psi*psi';
                psi = [randNum(5) + randNum(6)*1i; randNum(7)+randNum(8)*1i];
                psi = psi/norm(psi);
                B = psi*psi';
                psi = [randNum(9) + randNum(10)*1i; randNum(11)+randNum(12)*1i];
                psi = psi/norm(psi);
                C = psi*psi';
                pi = kron(A,kron(B,C));
                pi = reshape(pi,[64,1]);
                
                for j = 1:length(dt)
                    pi = U{j}*pi;
                    eigVal = eig(reshape(pi,[obj.dim^3,obj.dim^3]),'vector');
                    trRho(i,j) = sum(eigVal);
                    minEig(i,j) = min(real(eigVal));
                    maxImag(i,j) = max(imag(eigVal));
                end
            end
            trRho = sign(max(trRho-1)).*abs(max(trRho-1))+1;
            minEig = min(minEig);
            maxImag = max(maxImag);
        end
       
        function td = traceDistance(obj,state) %reduced state of spin 3
            td = 0.5*sum(abs(eig(obj.rho-state)));
        end
        
        function r = getYourMother(obj) %Doesn't make much sense...
            for i = 1:length(obj.coeffLOps)
                r(i) = obj.coeffLOps(i) * norm(obj.LOps);
            end
        end

        function r = getMaxDissRates(obj) 
            %gets the largest expected dissipation rates of all Lindblads,
            %using HS-norm, with sign if coeff coefficient is negative.
            for i = 1:length(obj.coeffLOps)
                %norm of operator gives greatest singular (or eigen) value
                r(i) = sign(obj.coeffLOps(i)) * norm( (obj.LOps{i})' * obj.LOps{i} );
            end
        end
        
        function r = getDissRates(obj)
            %gets the mean expected dissipation rates of all Lindblads,
            %using HS-norm, with sign if coeff coefficient is negative.
            for i = 1:length(obj.coeffLOps)
                %norm of operator gives greatest singular (or eigen) value
                r(i) = sign(obj.coeffLOps(i)) * trace( (obj.LOps{i})' * obj.LOps{i} );
            end
        end
        
        function p = checkState(obj) %get all eigenvalues and the trace
            p = eig(obj.rho);
            p = [p; sum(p)];
        end
        
        function getLindblad(obj, X, kappa, T) %Core routine to get the Lindblad operators for each bath
            %do this above at init
            %[eigvec, eigval] = eig(obj.H,'vector');
            X = sparse(X);
            
            expectX = obj.U'*X*obj.U;
            [r,c,s] = find(expectX);
            Egap = obj.E(r)-obj.E(c);
            [Em,~,idx] = uniquetol(Egap,10*eps);
            Em = Em.';
            K = cell(length(Em),1);
            for j = 1:length(Em)
                K{j} = sparse(zeros(obj.dim^3));
            end
            for k = 1:length(Egap)
                K{idx(k)} = K{idx(k)} + s(k) * obj.U(:,r(k))*obj.U(:,c(k))';
            end
           
%             Em = [];
%             K = cell(0);
%             for j = 1:length(obj.U)
%                 for k = 1:length(obj.U)
%                     expectX = obj.U(:,j)'*X*obj.U(:,k);
%                     if (abs(expectX)<10*eps)
%                         continue;
%                     end
%                     Egap = obj.E(j)-obj.E(k);
%                     idx = find(abs(Em-Egap)<eps*10);
%                     if (isempty(idx))
%                         Em = [Em Egap];
%                         K{end+1} = (obj.U(:,j)*obj.U(:,k)')*expectX;
%                     else
%                         K{idx} = K{idx}+(obj.U(:,j)*obj.U(:,k)')*expectX;
%                     end
%                 end
%             end
%             
            switch obj.meType
                case 'global'
                    for j = 1:length(Em)
                        W = Em(j);
                        switch obj.bathType
                            case 'ohmic'
                                gamma = abs(W)*(nbar(obj,abs(W),T) + heaviside(-W));
                            case 'flat'
                                gamma = (nbar(obj,abs(W),T) + heaviside(-W));
                            case 'cubic'
                                gamma = abs(W).^3*(nbar(obj,abs(W),T) + heaviside(-W));
                        end
                        obj.LOps{end+1} = sqrt(gamma*kappa)*K{j};
                    end
                    obj.bathIdx = [obj.bathIdx length(obj.LOps)];
                    return;
            end
%             gamma2 = zeros(length(Em));
%             for j = 1:length(Em)
%                 W = (Em(j) + Em )/2;
%                 absW = abs(W);
%                 idx = absW>10*eps;
%                 omega = (Em(j) - Em)*obj.deltat/2;
%                 gamma2(j,idx) =  absW(idx).*(nbar(obj,absW(idx),T) + heaviside(-W(idx))).*exp(1i*omega(idx)).*sinc(omega(idx)/pi);
%                 gamma2(j,~idx) = T * exp(1i*omega(~idx)).*sinc(omega(~idx)/pi);
%             end
            Emm = ones(length(Em),1)*Em;
            W = (Emm.' + Emm)/2;
            omega = (Emm.' - Emm)*obj.deltat/2;
            absW = abs(W);
            gamma = (absW.*(1./(exp(absW./T)-1) + heaviside(-W)));
%             gamma(absW<1e-3*T) = T + heaviside(-W(absW<1e-3*T)) ;
            gamma(absW<10*eps) = T;
            gamma = gamma.*exp(1i*omega).*sinc(omega/pi);

            
            %abs(gamma)
            [evec, ev] = eig(gamma,'vector');
            obj.eigVal = [obj.eigVal; ev];
            obj.coeffLOps = [obj.coeffLOps; ev];
%             %%%%%%%
%             %%%%%%%%%
            for j = 1:length(Em)
                obj.LOps{end+1} = zeros(obj.dim^3);
                for k = 1:length(Em)
                    obj.LOps{end} = obj.LOps{end} + evec(k,j)*K{k};
                end
                obj.LOps{end} = sqrt(kappa*ev(j))*obj.LOps{end};
            end
            obj.bathIdx = [obj.bathIdx length(obj.LOps)];
        end
             
    end
        
end


function [A, B] = ptrace(AB, dimA, dimB) %partial trace to get bipartition reduced states

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
