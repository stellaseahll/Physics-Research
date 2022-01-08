%2 COLLECTIVE spins coupled to two baths

classdef twoSpinWire < handle
    
    properties
        N1
        N2
        J1
        J2
        d %total dimension
        
        J1mat %{x,y,z,p,m,id} collective spin matrix cell size: 1x6
        J2mat
        
        T1 % unitless kB*T/hbar*omega
        T2
        n1
        n2
        
        w1 %should be finite for correct heat flow!
        w2
        k1
        k2
        g
        
        H %Hamiltonian with RESONANT exchange interaction
        L1 %Liouvillian with LOCAL thermal dissipator 1
        L2 %Liouvillian with LOCAL thermal dissipator 2
        %NB: L1' gives the adjunct Liouvillian applied to observables (I think...)
        L %total Liouvillian
        rhoSS %steady state (stored as Liouville vector)
    end
    
    methods
        function obj = twoSpinWire(Ns,G,ks,ns,ws)
            %init object.
            % Specify: 
            % 2-array of qubit numbers Ns = [N1,N2]
            % 1 coupling rate g
            % thermalization rates ks = [k1,k2]
            % average occupation numbers ns = [n1,n2]
            % (optional) eigenfrequencies ws = [w1,w2] (default:[1,1], i.e.
            %   resonant model. WARNING! w1,w2 != 0 for heat flows!
            obj.N1 = Ns(1); 
            obj.N2 = Ns(2);
            obj.g = G;
            obj.k1 = ks(1);
            obj.k2 = ks(2);
            obj.n1 = ns(1);
            obj.n2 = ns(2);
            if nargin==5
                obj.w1 = ws(1);
                obj.w2 = ws(2);
            else
                obj.w1 = 1;
                obj.w2 = 1;
            end
            obj.initOperators();
        end
        
        function initOperators(obj)
            %init all operators and derived quantities
            obj.J1 = obj.N1/2;
            obj.J2 = obj.N2/2;
            obj.d = (obj.N1+1)*(obj.N2+1);
            obj.T1 = 1/log(1+1/obj.n1);
            obj.T2 = 1/log(1+1/obj.n2);
            
            [jz,jp,jm,jx,jy] = collSpinOps(obj.J1);
            obj.J1mat = {jx,jy,jz,jp,jm,speye(obj.N1+1)};
            [jz,jp,jm,jx,jy] = collSpinOps(obj.J2);
            obj.J2mat = {jx,jy,jz,jp,jm,speye(obj.N2+1)};
            
            obj.H = obj.w1*kron(obj.J1mat{3},obj.J2mat{6}) + obj.w2*kron(obj.J1mat{6},obj.J2mat{3}) ...
                + obj.g*( kron(obj.J1mat{4},obj.J2mat{5}) + kron(obj.J1mat{5},obj.J2mat{4}) );
            
            obj.L1 = obj.k1*(obj.n1+1)*spDissipator( kron(obj.J1mat{5},obj.J2mat{6}) ) ...
                + obj.k1*obj.n1*spDissipator( kron(obj.J1mat{4},obj.J2mat{6}) );
            obj.L2 = obj.k2*(obj.n2+1)*spDissipator( kron(obj.J1mat{6},obj.J2mat{5}) ) ...
                + obj.k2*obj.n2*spDissipator( kron(obj.J1mat{6},obj.J2mat{4}) );
            
            obj.L = 1i * ( spLeftMultiply(obj.H) - spRightMultiply(obj.H) ) ...
                + obj.L1 + obj.L2;
        end
        
        function steadyState(obj)
            obj.rhoSS = spnull(obj.L);
            ss = size(obj.rhoSS);
            if ss(2)>1
                warning(' %i null vectors found!',ss(2));
            end
            %obj.rhoSS = reshape(obj.rhoSS(:,1),[obj.d,obj.d]);
            %normalize!
            trRho = sum(obj.rhoSS(1:obj.d+1:end,:));
            obj.rhoSS = obj.rhoSS ./ trRho;
        end
        
        function [Q1,Q2] = getHeatFlows(obj)
            %computes the heat flows from bath 1 and 2 at steady state
            h = obj.H(:).'; %vectorized observable
            Q1 = full( h * (obj.L1*obj.rhoSS) );
            Q2 = full( h * (obj.L2*obj.rhoSS) );
        end
        
    end
        
end