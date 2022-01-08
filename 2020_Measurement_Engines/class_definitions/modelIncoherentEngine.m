% 10/5/2019 Incoherent feedback engine with qubit + HO pointer
% Projective meas. in displaced Fock basis. Use only up to Fock number
% where potential energy crosses zero (but at least ground state)
% Basis: |g><g||R_m><R_m|, |e><e||L_m><L_m|
% |R_m> = D(g/2wp) |m>, |L_m> = D(-g/2wp) |m>, |m> = Fock
classdef modelIncoherentEngine < handle    
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
        % Interaction strengths
        km % measurement rate
        g % interaction between spin and pointer
        
        %Hamiltonians
        Hq % spin
        Hp % pointer
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
        spLm %measurement 
        spL %total
        
        %Dissipators
        Lh % hot bath
        Lc % cold bath
        Lm %measurement (contains projector w/ and w/o bit flip)
        
        %Other operators
        b % b = a + g/wp sz, diagonal in basis used
        C %\sum_m,n |R_m><L_n| <R_m|L_n>
        %CAREFUL: Product operators qubit-osci are problematic in
        %convergence when using the conditionally displaced basis. Products
        %must be formulated separately including the correct overlap matrix
        %elements of oscillator states (obj.C), otherwise convergence
        %issues. Reason: Poor convergence of overlap integrals of displaced
        %Fock states with cutoff.
        
    end
    
    methods
        
        function obj = modelIncoherentEngine(dimp,ws,wp,g,km,kh,kc,Th,Tc,dissType)
            obj.dimp = dimp;
            obj.ws = ws;
            obj.wp = wp;
            obj.g = g;
            obj.km = km;
            obj.kh = kh;
            obj.kc = kc;
            obj.Th = Th;
            obj.Tc = Tc;

            obj.createOperators();
            obj.createHamiltonians();
            obj.createDissipators(dissType);
            obj.createLiouville();
        end

        function createOperators(obj) 
            obj.createb();
            obj.createC();
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

        function createDissipators(obj,dissType) %get Lindblad operators            
            obj.createHotDissipators(dissType);
            obj.createColdDissipators();
            obj.createMeasurementP();
        end
        
        function createColdDissipators(obj) %get Lindblad operators
            obj.nc = 1/(exp(obj.wp/obj.Tc)-1); 
            obj.Lc{1} = sqrt(obj.kc*(obj.nc+1))*obj.b;
            obj.Lc{2} = sqrt(obj.kc*obj.nc)*obj.b';
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
            %global dissipators: incoherent in displaced Fock basis
            %transitions
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
        
        function createMeasurementP(obj) %creates the projector
            %compute max Fock number of projector based on potential energy
            %at x=0 assuming x=(a+a')/sqrt(2) i.e. Hp=wp*(P^2/2+X^2/2)
            Nmax = max( [(obj.g/obj.wp)^2/4 - 0.5, 1] ); %not integer
            %compute the projector
            u = ones(obj.dimp,1);
            u( (ceil(Nmax)+1):end )=0;
            %projector w/o bit flip
            P = kron([0,0;0,1],diag(u));
            %P = P + kron([1,0;0,0], obj.C*(diag(u)*obj.C'));
            %sxP = kron([0,0;1,0],diag(u)*obj.C');
            P = P + kron([1,0;0,0], obj.C'*(diag(u)*obj.C));
            sxP = kron([0,0;1,0],diag(u)*obj.C);
            
            sxP = sxP + sxP';
            obj.Lm{1} = sqrt(obj.km)*sxP;
            obj.Lm{2} = sqrt(obj.km)*P; %is id-P, but doesn't matter for dissipator
        end
        
        function createLiouville(obj) %Liouville superoperator from unitary and Lindblads
            obj.createUnitary();
            obj.spLh = obj.createL(obj.Lh);
            obj.spLc = obj.createL(obj.Lc);
            obj.spLm = obj.createL(obj.Lm);
            obj.spL = obj.spU + obj.spLh + obj.spLc + obj.spLm;
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
       
        function updateColdTemp(obj,T)
            obj.Tc = T;
            obj.createColdDissipators();
            obj.spLc = obj.createL(obj.Lc);
            obj.spL = obj.spU + obj.spLh + obj.spLc;
        end
        
        function updateHotTemp(obj,T)
            obj.Th = T;
            obj.createHotDissipators(dissType);
            obj.spLh = obj.createL(obj.Lh);
            obj.spL = obj.spU + obj.spLh + obj.spLc;
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
            obj.H = obj.Hq + obj.Hp;
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
        
        function [Qm, Qh, Qc, Qba] = getThermoProp(obj)
            %careful: Qm contains work and backaction! It is total energy
            %change due to measurement dissipator
            Qm = obj.getQ(obj.Lm,obj.H);
            Qh = obj.getQ(obj.Lh,obj.H);
            Qc = obj.getQ(obj.Lc,obj.H);
            %measurement backaction: total energy change due to measurement
            %w/o feedback! NB: Lm{2} is D[P]
            Qba = 2*trace(obj.diss(obj.sprhoSS,obj.Lm{2})*obj.H);
        end
        
        function [Qm, Qh, Qc] = getThermoPropBare(obj)
            %Computes only change of bare qubit energy!
            Qm = obj.getQ(obj.Lm,obj.Hq);
            Qh = obj.getQ(obj.Lh,obj.Hq);
            Qc = obj.getQ(obj.Lc,obj.Hq);
        end
        
        function Q = getQ(obj,L,H)
            Q = 0;
            for i = 1:length(L)
                Q = Q + trace(obj.diss(obj.sprhoSS,L{i})*H);
            end
        end
        
        function d = diss(obj,rho,L)
            LdgL = L'*L;
            d = L*rho*L' - 0.5*(rho*LdgL + LdgL*rho);
        end
% 
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
