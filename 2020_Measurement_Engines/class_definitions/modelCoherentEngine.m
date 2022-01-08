% 9/5/2019 Coherent feedback engine with qubit + HO pointer
% Basis: |g><g||R_m><R_m|, |e><e||L_m><L_m|
% |R_m> = D(g/2wp) |m>, |L_m> = D(-g/2wp) |m>, |m> = Fock
classdef modelCoherentEngine < handle    
    properties
        % Dimension of pointer
        dimp
        % Temperature of bath
        Th
        Tc
        % Bath occupation
        nh %vector
        nc
        % Bath coupling
        kh
        kc
        % Frequencies
        wp % pointer
        ws % spin
        wd % driving field
        % Interaction strengths
        G % feedback
        g % interaction between spin and pointer
        
        %Hamiltonians
        Hq % spin
        Hp % pointer
        Hd % driving
        Hrot % rotated frame
        H %total
        
        %States
        sprhoSS %total
        sprhoQ %reduced spin
        sprhoR %reduced right pointer state
        sprhoL %reduced left pointer state
        
        
        %Superoperators
        spU %unitary
        spLh %hot bath 
        spLc %cold bath 
        spL %total
        
        %Dissipators
        Lh % hot bath
        Lc % cold bath
        
        %Other operators
        b % b = a + g/wp sz, diagonal in basis used
        C %\sum_m,n |R_m><L_n| <R_m|L_n>
        F % Use if there is a Gaussian fx coupling accompanying feedback Ham
        %CAREFUL: Product operators qubit-osci are problematic in
        %convergence when using the conditionally displaced basis. Products
        %must be formulated separately including the correct overlap matrix
        %elements of oscillator states (obj.C), otherwise convergence
        %issues. Reason: Poor convergence of overlap integrals of displaced
        %Fock states with cutoff.
        
        dissType
    end
    
    methods
        
        function obj = modelCoherentEngine(dimp,ws,wp,g,G,kh,kc,Th,Tc,f,dissType,delta)
            obj.dimp = dimp;
            obj.ws = ws;
            obj.wp = wp;
            obj.g = g;
            obj.G = G;
            obj.kh = kh;
            obj.kc = kc;
            obj.Th = Th;
            obj.Tc = Tc;
            obj.dissType = dissType;
            if (nargin == 12)
                obj.wd = ws-delta;
            else
                obj.wd = ws - g*g/wp; %if no detuning set, assume optimal of delta = g*g/wp
            end
            obj.createOperators(f);
            obj.createHamiltonians();
            obj.createDissipators(dissType);
            obj.createLiouville();
        end

        function createOperators(obj,f)            
            obj.createb();
            obj.createC();
            obj.createF(f)
        end
        
        function createb(obj)
            %in the conditionally displaced Fock basis, the displaced mode
            %operator is just the undisplaced one.
            a = diag(sqrt(1:obj.dimp-1),1);
            obj.b = sparse(kron(eye(2),a));
        end
        
        function createC(obj)
            %Creates matrix C_{mn} w/ elements <m|D(g/wp)|n> = <Rm|Ln>
            a = diag(sqrt(1:(obj.dimp-1)),1); 
            %direct calculation, guarantees unitarity.
            obj.C = expm(obj.g/obj.wp*(a'-a));
%             obj.C = zeros(obj.dimp);
%             for i = 0:obj.dimp-1
%                 for j = i:obj.dimp-1
%                     obj.C(i+1,j+1) = obj.getmDn(j,i,obj.g/obj.wp);
%                     obj.C(j+1,i+1) = obj.C(i+1,j+1);
%                 end
%             end
        end

        function k = getmDn(obj,m,n,alpha) 
            %OBSOLETE: Has some weird sign problem, which leads to nonunitary operator
            k = 0;
            for l = 0:n
                k = k + (-1)^l *nchoosek(m,n-l)*alpha^(2*l)/gamma(l+1);
            end
            k = k*sqrt(gamma(n+1)/gamma(m+1))*alpha^(m-n)*exp(-alpha^2/2);
        end

        function createF(obj,f)
            if ~isa(f, 'function_handle') 
                obj.F = obj.C;
                return;
            end
            obj.F = zeros(obj.dimp);
            dx = obj.g/sqrt(2)/obj.wp;
            %new implementation of H.O. eigenstates
            Hcoeff = myHermiteCoeff(obj.dimp-1); %coefficient matrix for Hermite poly's
            %eigenstates w/o pi^{-1/4}
            psiHO = @(nn,x) polyval( Hcoeff(nn+1,1:(nn+1)), x ).*exp(-x.*x/2);
            for m = 0:(obj.dimp-1)
                for n = 0:(obj.dimp-1)
                     %obj.F(m+1,n+1) = integral(@(x) obj.myHermite(m,x+dx).*obj.myHermite(n,x-dx).*f(x).*exp(-(x+dx).^2/2).*exp(-(x-dx).^2/2),-inf,inf)/sqrt(2^(m+n)*factorial(m)*factorial(n)*pi);    
                     %numerical guess of reasonable integral bounds
                     lims = (dx + 2*max([ sqrt(2*n+1), sqrt(2*m+1), 5 ]));
                     %lims=inf;
                     obj.F(m+1,n+1) = integral(@(x) psiHO(m,x+dx).*psiHO(n,x-dx).*f(x),-lims,lims)/sqrt(pi);    
                end
            end
        end
        
        function createDissipators(obj,dissType) %get Lindblad operators           
            obj.createHotDissipators(dissType);
            obj.createColdDissipators();
        end
        
        function createColdDissipators(obj) %get Lindblad operators
            obj.nc = 1/(exp(obj.wp/obj.Tc)-1); 
            obj.Lc{1} = sqrt(obj.kc*(obj.nc+1))*obj.b;
            obj.Lc{2} = sqrt(obj.kc*obj.nc)*obj.b';
        end
        
        function updateColdTemp(obj,T)
            obj.Tc = T;
            obj.updateColdDissipators;
        end
        
        function updateHotTemp(obj,T)
            obj.Th = T;
            obj.updateHotDissipators;
        end
        
        function updateColdThermRate(obj,kc)
            obj.kc = kc;
            obj.updateColdDissipators;
        end
        
        function updateHotThermRate(obj,kh)
            obj.kh = kh;
            obj.updateHotDissipators;
        end
        
        function updateHotDissipators(obj)
            obj.createHotDissipators(obj.dissType);
            obj.spLh = obj.createL(obj.Lh);
        end
        
        function updateColdDissipators(obj)
            obj.createColdDissipators();
            obj.spLc = obj.createL(obj.Lc);
        end
        
        function updateLiouville(obj)
            obj.spL = obj.spU + obj.spLh + obj.spLc;
        end
        
        function updateIntStrength(obj,g)
            obj.g = g;
            obj.createC();
            obj.createF();
            obj.createHamiltonians();
            obj.createUnitary();
        end
        
        
        function createHotDissipators(obj,dissType) %get Lindblad operators
            if (dissType == 'l')
                obj.nh = 1/(exp(obj.ws/obj.Th)-1);
                %obj.Lh{1} = sqrt(obj.kh*(obj.nh+1))*kron([0 1; 0 0],obj.C);
                %obj.Lh{2} = sqrt(obj.kh*obj.nh)*kron([0 0; 1 0],obj.C');
                % |e><g| matrix elements are that of D^2 = C
                obj.Lh{1} = sqrt(obj.kh*(obj.nh+1))*kron([0 1; 0 0],obj.C');
                obj.Lh{2} = sqrt(obj.kh*obj.nh)*kron([0 0; 1 0],obj.C);
                return;
            end
            %global: Not strictly global, since driving Hamiltonian is
            %omitted! Invalid for strong driving...
            X = cell(2*obj.dimp+1);
            E = zeros(1,2*obj.dimp+1);
            X{1} = sparse(kron([0 0; 1 0],diag(diag(obj.C))));
            E(1) = obj.ws;
            for k = 1:obj.dimp
                E(1+k) = obj.ws+obj.wp*k;
                X{1+k} = sparse(kron([0 0; 1 0],diag(diag(obj.C,-k),-k)));
                E(1+k+obj.dimp) = obj.ws-obj.wp*k;
                X{1+k+obj.dimp} = sparse(kron([0 0; 1 0],diag(diag(obj.C,k),k)));
            end
            
            obj.nh = 1./(exp(E./obj.Th)-1);
            for i = 1:length(X)
                obj.Lh{i} = sqrt(obj.kh*(obj.nh(i)+1))*X{i}';
                obj.Lh{length(X)+i} = sqrt(obj.kh*obj.nh(i))*X{i};
            end
        end
        
        function createLiouville(obj) %Liouville superoperator from unitary and Lindblads
            obj.createUnitary();
            obj.spLh = obj.createL(obj.Lh);
            obj.spLc = obj.createL(obj.Lc);
            obj.spL = obj.spU + obj.spLh + obj.spLc;
        end
        
        function createUnitary(obj) 
            obj.spU = -1i* ( obj.spLeftMultiply(obj.H) - obj.spRightMultiply(obj.H) );
        end
        
        function L = createL(obj,LOps)
            L = sparse(4*obj.dimp^2,4*obj.dimp^2);
            for i = 1:length(LOps)
                L = L + obj.computeDissSuperOp(LOps{i});
            end
        end
       
        function f=myHermite(obj,m,x)
            if (m==0)
                f = 1;
                return;
            end
            if (m==1)
                f = 2*x;
                return;
            end
            Hn(1,1:length(x)) = 1;
            Hn(2,1:length(x)) = 2*x;
            for i = 3:(m+1)
                Hn(i,1:length(x)) = 2*x.*Hn(i-1,1:length(x))-2*(i-2)*Hn(i-2,1:length(x));
            end
            f = Hn(end,1:length(x));
        end
        
        function createHamiltonians(obj)
            sz = sparse(kron([-1 0; 0 1],eye(obj.dimp)));
            obj.Hq = obj.ws/2*sz;
            obj.Hp = obj.wp*obj.b'*obj.b;
            obj.Hd = sparse(obj.G*(kron([0 0;1 0],obj.F) +kron([0 1; 0 0],obj.F')));
            obj.Hrot = obj.wd/2*sz;
            obj.H = obj.Hq + obj.Hp + obj.Hd - obj.Hrot;
        end
        
        function M = computeDissSuperOp(obj,A)
            if ~issparse(A)
                A = sparse(A);
            end
            AdgA = A'*A;
            M = obj.spLrMultiply(A) - 0.5*obj.spLeftMultiply(AdgA) - 0.5*obj.spRightMultiply(AdgA);
        end
        
        function M = spLrMultiply(obj,A)
            M = (kron(conj(A),A));
        end

        function M = spLeftMultiply(obj,A)
            M = kron(speye(size(A)),A); %if one of the matrices is sparse, kron spits out sparse
        end
        
        function M = spRightMultiply(obj,A)
            M = kron(A.',speye(size(A))); %if one of the matrices is sparse, kron spits out sparse
        end
        
        function findSS(obj) %finds the steady state of L
            obj.sprhoSS = spnull(obj.spL);
            obj.sprhoSS = reshape(obj.sprhoSS,obj.dimp*2,obj.dimp*2);
            obj.sprhoSS = obj.sprhoSS/trace(obj.sprhoSS);
            obj.sprhoR = obj.sprhoSS(1:obj.dimp,1:obj.dimp);
            obj.sprhoL = obj.sprhoSS((1:obj.dimp)+obj.dimp,(1:obj.dimp)+obj.dimp);
            obj.sprhoQ = [trace(obj.sprhoR) 0; 0 trace(obj.sprhoL)];
        end
        
        function [W Qh Qc] = getThermoProp(obj)
            W =obj.getWork();
            Qh = obj.getQ(obj.Lh,obj.H + obj.Hrot);
            Qc = obj.getQ(obj.Lc,obj.H + obj.Hrot);
        end
        
        function W = getWork(obj)
           W = obj.G*obj.wd*trace((-kron([0 0; 1i 0],obj.F) + kron([0 1i;0  0],obj.F'))*obj.sprhoSS);
           W = real(W);
        end
        
        function Q = getQ(obj,L,H)
            Q = 0;
            for i = 1:length(L)
                Q = Q + trace(obj.diss(obj.sprhoSS,L{i})*H);
            end
            Q = real(Q);
        end
        
        function d = diss(obj,rho,L)
            LdgL = L'*L;
            d = L*rho*L' - 0.5*(rho*LdgL + LdgL*rho);
        end

        
%         function [A, B] = ptrace(AB, dimA, dimB)
% 
%             B = zeros(dimB);
%             A = zeros(dimA);
%             for i = 0:dimA-1
%                 for j = 0:dimA-1
%                     if (i==j)
%                         tmp = AB(i*dimB+1:(i+1)*dimB,i*dimB+1:(i+1)*dimB);
%                         A(i+1,i+1) = trace(tmp);
%                         B = B + tmp;
%                     else
%                        tmp = AB(i*dimB+1:(i+1)*dimB,j*dimB+1:(j+1)*dimB);
%                        A(i+1,j+1) = trace(tmp);
%                     end
%                 end 
%             end
%         end


       
            
    end
        
end
