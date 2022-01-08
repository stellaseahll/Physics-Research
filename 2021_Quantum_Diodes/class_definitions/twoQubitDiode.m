
classdef twoQubitDiode < handle
    properties
        %eigenfrequencies of qubits
        E1=1;
        E2=2;
        %coupling between qubits
        g=0.01; %3-spin coupling rate (only resonant RWA terms!)
        %bath temperature
        T1 = 2.0;
        T2 = 1.0;
        %bath coupling
        k1 = 1e-3;
        k2 = 1e-3;
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
        eigVal = []; %eigenvalues of gamma matrix
        
        %bath type
        bathType = 'flat'; %bath type: 'flat' or 'ohmic'
        intType = 'RWA';  %interaction type: 'RWA' or 'full'
        meType = 'partial'; %master equation type: 'partial','local','global','breuer'
        %states
        rho % steady state, i.e. kernel element of L
        rho1
        rho2
        %permutation
        %coarsegrain time
        deltat
        gamma = cell(1,3)
        
        %% Steady heat current
        Q
        %% standard operators
        sz1 = sparse(kron([-0.5 0; 0 0.5],eye(2)));
        sz2 = sparse(kron(eye(2),[-0.5 0; 0 0.5]));
        sx1 = sparse(kron([0 1; 1 0],eye(2)));
        sx2 = sparse(kron(eye(2),[0 1; 1 0]));
        sp1 = sparse(kron([0 0; 1 0],eye(2)));
        sp2 = sparse(kron(eye(2),[0 0; 1 0]));
        sm1 = sparse(kron([0 1; 0 0],eye(2)));
        sm2 = sparse(kron(eye(2),[0 1; 0 0]));
        hint = sparse([0 0 0 0; 0 0 1 0;0 1 0 0; 0 0 0 0]);
    end
    
    methods
        
        function obj = twoQubitDiode(g,T,k,E,Btype,Mtype)
            if nargin>=1
                obj.g = g;
            end
            if nargin>=2
                obj.T1 = T(1);
                obj.T2 = T(2);
            end
            if nargin>=3
                obj.k1 = k(1);
                obj.k2 = k(2);
            end
            if nargin>=4
                obj.E1 = E(1);
                obj.E2 = E(2);
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
            obj.deltat = 1./sqrt(max([obj.E1,obj.E2]));
            obj.resetOperators();
            obj.findSS();
            obj.getTotalHeatFlow();
        end
        %
        %         function tmp = getOp(obj,eigvec,i,j,X)
        %             tmp = (eigvec(:,i)*eigvec(:,j)') * (eigvec(:,i)'*X*eigvec(:,j));
        %         end
        
        
        
        function resetOperators(obj) %create all operators and superops.
            obj.Hint = obj.g*obj.hint;
            obj.H = obj.E1*obj.sz1 + obj.E2*obj.sz2 + obj.Hint;
            [obj.U , obj.E] = eig(full(obj.H),'vector');
            obj.U(abs(obj.U)<10*eps) = 0;
            obj.U = sparse(obj.U);
            obj.L = -1i*obj.spLeftMultiply(obj.H) + 1i*obj.spRightMultiply(obj.H);
            switch obj.meType
                case 'local'
                    obj.bathIdx = [0 2 4];
                    n1 = nbar(obj,obj.E1,obj.T1);
                    n2 = nbar(obj,obj.E2,obj.T2);
                    switch obj.bathType
                        case 'cubic'
                            K1 = obj.k1*(obj.E1)^3;
                            K2 = obj.k2*(obj.E2)^3;
                        case 'ohmic'
                            K1 = obj.k1*(obj.E1);
                            K2 = obj.k2*(obj.E2);
                        case 'flat'
                            K1 = obj.k1;
                            K2 = obj.k2;
                    end
                    obj.LOps{1} = sqrt(K1*(n1+1)) * obj.sm1;
                    obj.LOps{2} = sqrt(K1*n1) * obj.sp1;
                    obj.LOps{3} = sqrt(K2*(n2+1)) * obj.sm2;
                    obj.LOps{4} = sqrt(K2*n2) * obj.sp2;
                otherwise
                    obj.bathIdx = 0;
                    getLindblad(obj,obj.sx1,obj.k1,obj.T1);
                    getLindblad(obj,obj.sx2,obj.k2,obj.T2);
            end
            %add the Lindblad terms
            for j = 2:length(obj.bathIdx)
                a=obj.bathIdx(j-1)+1;
                b=obj.bathIdx(j);
                for k = a:b
                    LL = sparse(obj.LOps{k});
                    tmp = LL'*LL;
                    obj.L = obj.L + obj.spLrMultiply(LL) - ...
                        0.5*( obj.spLeftMultiply(tmp) + obj.spRightMultiply(tmp));
                end
            end
        end

        
        function SS = findSS(obj) %finds the steady state of L
            SS = obj.spnull(obj.L);
            SS = reshape(full(SS),[4,4]);
            SS = SS/trace(SS);
            obj.rho = SS;
            [obj.rho1, obj.rho2] = obj.ptrace(SS,2,2);
        end
        
        function [rho1 rho2] = generateRhoTh(obj) %generate product of thermal states
            rho1 = exp(-obj.E1/obj.T1*(0:1));
            rho1 = diag(rho1/sum(rho1));
            rho2 = exp(-obj.E2/obj.T2*(0:1));
            rho2 = diag(rho2/sum(rho2));
        end
        
        
        function U=generateU(obj,dt)
            U = expm(obj.L*dt);
        end
        
        %% MATH OPERATIONS
        function M = spLeftMultiply(obj,A)
            M = kron(speye(size(A)),A); %if one of the matrices is sparse, kron spits out sparse
        end
        
        function M = spRightMultiply(obj,A)
            M = kron(A.',speye(size(A))); %if one of the matrices is sparse, kron spits out sparse
        end
        
        function M = spLrMultiply(obj,A)
            if ~issparse(A)
                A = sparse(A);
            end
            M = (kron(conj(A),A));
        end
        
        function Z = spnull(obj, S, varargin)
            if issparse(S)
                [m n] = size(S);
                try
                    %[Q R E] = qr(S.'); %#ok %full QR
                    %SN 23.04.19: For complex matrices, one should probably use dagger
                    %instead of transpose
                    [Q R E] = qr(S'); %#ok %full QR
                    if m > 1
                        s = diag(R);
                    elseif m == 1
                        s = R(1);
                    else
                        s = 0;
                    end
                    s = abs(s);
                    tol = norm(S,'fro') * eps(class(S));
                    r = sum(s > tol);
                    Z = Q(:,r+1:n);
                catch %#ok
                    % sparse QR is not available on old Matlab versions
                    err = lasterror(); %#ok
                    if strcmp(err.identifier, 'MATLAB:maxlhs')
                        Z = null(full(S), varargin{:});
                    else
                        rethrow(err);
                    end
                end
            else % Full matrix
                Z = null(S, varargin{:});
            end
        end
        
        
        function [A, B] = ptrace(obj,AB, dimA, dimB)
            
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
        
        
        
        %% THERMO QUANTITIES
        
        function n = nbar(obj,omega,T)
            n = 1./(exp(omega./T)-1);
            %cutoff criterion
            %             n(T>100*omega) = T./omega(T>100*omega);
        end
        
        
        function Q = getHeatFlow(obj,j) %average energy change at SS coming from Lindblad term j
            %CHECK, modify as needed
            LL = obj.LOps{j};
            Q = trace( obj.H * ( LL*(obj.rho*(LL')) - 0.5*(LL'*LL)*obj.rho - 0.5*obj.rho*(LL'*LL) ) );
            
        end
        
        
        function getTotalHeatFlow(obj)
            for bath = 1:2
                a = obj.bathIdx(bath)+1;
                b = obj.bathIdx(bath+1);
                obj.Q(bath) = 0;
                for k = a:b
                    obj.Q(bath) = obj.Q(bath) + getHeatFlow(obj,k);
                end
                obj.Q(bath) = real(obj.Q(bath));
            end
        end
        
        
        
        function E = getEnergy(obj,bath) %avg single-mode energy of spin 'bath' in SS
            switch bath
                case 1
                    E = obj.rho1(2,2)*obj.E1;
                case 2
                    E = obj.rho2(2,2)*obj.E2;
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
                    eigVal = eig(reshape(pi,[4,4]),'vector');
                    trRho(i,j) = sum(eigVal);
                    minEig(i,j) = min(real(eigVal));
                    maxImag(i,j) = max(imag(eigVal));
                end
            end
            trRho = sign(max(trRho-1)).*abs(max(trRho-1))+1;
            minEig = min(minEig);
            maxImag = max(maxImag);
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
                K{j} = sparse(zeros(4));
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
                obj.LOps{end+1} = zeros(4);
                for k = 1:length(Em)
                    obj.LOps{end} = obj.LOps{end} + evec(k,j)*K{k};
                end
                obj.LOps{end} = sqrt(kappa*ev(j))*obj.LOps{end};
            end
            obj.bathIdx = [obj.bathIdx length(obj.LOps)];
        end
        
    end
    
end