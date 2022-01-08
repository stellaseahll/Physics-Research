% Class definition for N qubits interacting with N
% qubits.

classdef spEqualSpins < spCollectiveSpins
    
    properties
    end
    
    methods
        function obj = spEqualSpins(N,T,k,w)
            if (nargin<=2)
                k = 0.001*[1 1];
            end
            if (nargin<=3)
                w = [1 1];
            end
            obj = obj@spCollectiveSpins([N N],T,k,w);
            
        end
        
        function addInteraction(obj,g,ops)
            %Input ops is a rx2 matrix that specifies the operators for the
            %interaction where 1=x,2=y,3=z,4=p,5=m,6=pm,7=mp
            
            [r, c] = size(ops);
            if (c~=2)
                error('Check input matrix for interaction.');
            end
            for k = 1:r
                obj.H = obj.H + g*(obj.J1mat{ops(k,1)}*obj.J2mat{ops(k,2)});
                for j = 1:obj.N1
                    obj.indH = obj.indH + g*(obj.indJ1mat{j,ops(k,1)}*obj.indJ2mat{j,ops(k,2)});
                end
            end
        end
        
        function updateInteraction(obj,g,ops)
            if (nargin==1)
                obj.getIndividualLiouville();
                obj.getCollectiveLiouville();
                return;
            end
            obj.addInteraction(g,ops);
            obj.getIndividualLiouville();
            obj.getCollectiveLiouville();
            obj.getAllProjectors();
        end
        
        function rhof = swapIndex(obj,rho,j,k,N)
            %swaps index j and k for N qubits
            j00 = kron(speye(2^(j-1)),kron(sparse([1 0; 0 0]),speye(2^(N-j))));
            j01 = kron(speye(2^(j-1)),kron(sparse([0 1; 0 0]),speye(2^(N-j))));
            j11 = kron(speye(2^(j-1)),kron(sparse([0 0; 0 1]),speye(2^(N-j))));
            k00 = kron(speye(2^(k-1)),kron(sparse([1 0; 0 0]),speye(2^(N-k))));
            k10 = kron(speye(2^(k-1)),kron(sparse([0 0; 1 0]),speye(2^(N-k))));
            k11 = kron(speye(2^(k-1)),kron(sparse([0 0; 0 1]),speye(2^(N-k))));
            tmp = j01*k10 ;
            U = speye(2^N) - j00*k11 - j11*k00 + tmp + tmp';
            rhof = U*rho*U';
        end
        
        function prepareNPairs(obj,rho)
            %prepares N pairs of qubits given by rho
            obj.rho0 = obj.createNPairs(rho,obj.N1);
        end
        
        function prepareAllIdentical(obj,rho)
            %prepares all qubits in S1 and S2 in rho
            obj.rho0 = obj.createNIdentical(rho,obj.N1*2);
        end
        
        function prepareSepIdentical(obj,rho1,rho2)
            %prepares all qubits in S1 in rho1^n and all qubits in S2 in
            %rho2^n
            obj.rho0 = kron(obj.createNIdentical(rho1,obj.N1),obj.createNIdentical(rho2,obj.N1));
        end
        
        function prepareSep(obj,rho1,rho2)
            %prepares all qubits in S1 in rho1 and all qubits in S2 in rho2
            [r1,c1] = size(rho1);
            [r2,c2] = size(rho2);
            if (any([r1 c1 r2 c2]~=2^obj.N1))
                error('Input state should be a 2-qubit %dx%d matrix.',2^obj.N1,2^obj.N1);
            end
            obj.rho0 = kron(rho1,rho2);
        end
        
        function Rho = createNPairs(obj,rho,N)
            %createNpairs creates N copies of 2 qubit state given by rho
            [r,c] = size(rho);
            if (r~=4&&c~=4)
                error('Input state should be a 2-qubit 4x4 matrix.');
            end
            rho = sparse(rho);
           
            Rho = rho;
            
            for i = 2:N
                Rho = kron(Rho,rho);
                idxStart=i;
                idxEnd=2*i-1;
                for j=idxEnd:-1:idxStart+1
                    Rho = obj.swapIndex(Rho,j,j-1,2*i);
                end
            end
        end
        
        function Rho = createNIdentical(obj,rho,N)
            %createNidentical creates N copies of a single qubit state
            %given by rho
            [r,c] = size(rho);
            if (r~=2&&c~=2)
                error('Input state should be a qubit 2x2 matrix.');
            end
            rho = sparse(rho);
            Rho = rho;
            for i = 2:N
                Rho = kron(Rho,rho);
            end
        end
        

    end
end