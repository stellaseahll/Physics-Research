%Class definition for repeated interaction at Poisson times, system and
%multiple baths. System is 1 spin, 1 bath

classdef collisionModelSpinBath < collisionModelThermal
    
    properties 
        %additional to parent class
        omegaS
        %spin operators
        JS
        %JzS
        %JxS
        %JyS
        %JpS
        %JmS
        JB
        %JzB
        %JxB
        %JyB
        %JpB
        %JmB
        %coupling coefficients
        %gxx = 0;
        %gyy = 0;
        %gzz = 0;
        %gpm = 0;
        gsb %matrix of coupling rates
        Lpm
    end
    
    methods
        function obj = collisionModelSpinBath(ds,omegas,db,omegab,Tb,gamma,dtint,gsb,rho0)
            %call constructor of parent class
            obj = obj@collisionModelThermal(ds,1,db,omegab,Tb,gamma,dtint);
            %do more specific stuff
            J = (ds-1)/2;
            m = (-J:J).';
            obj.omegaS = omegas;
            obj.ES = omegas*m;
            obj.gsb = gsb;
            %Assume we start in ground state
            obj.rhoS0 = rho0;
            %spin operators S
            %obj.JzS = diag(m);
            t = sqrt((J-m).*(J+m+1));
            JpS = diag(t(1:end-1),-1);
            %obj.JxS = 0.5*(obj.JpS + obj.JpS');
            %obj.JyS = 0.5*1i*(obj.JpS - obj.JpS'); %CHECK!
            obj.JS{1} = 0.5*(JpS + JpS'); %Jx
            obj.JS{2} = 0.5*1i*(JpS - JpS'); %Jy
            obj.JS{3} = diag(m); %Jz
            obj.JS{4} = JpS; %Jp
            obj.JS{5} = JpS'; %Jm
            %spin operators B
            J = (db-1)/2;
            m = (-J:J).';
            %obj.JzB = diag(m);
            t = sqrt((J-m).*(J+m+1));
            JpB = diag(t(1:end-1),-1);
            %obj.JxB = 0.5*(obj.JpB + obj.JpB');
            %obj.JyB = 0.5*1i*(obj.JpB - obj.JpB');
            obj.JB{1} = 0.5*(JpB + JpB'); %Jx
            obj.JB{2} = 0.5*1i*(JpB - JpB'); %Jy
            obj.JB{3} = diag(m); %Jz
            obj.JB{4} = JpB; %Jp
            obj.JB{5} = JpB'; %Jm
            %prepare interaction Hamiltonian
            obj.prepareHint(gsb);
% 			obj.getLshort([gsb gsb(end)]*dtint,gamma);
            
%             obj.getInsult();
        end
        
%         function prepareHint(obj,g) 
%             %must be called before prepareSim of parent class
%             %Hint must be represented in eigenbasis of HS!!!
%             obj.Hint{1} = g(4)*kron(obj.JpB,obj.JpS');
%             obj.Hint{1} = obj.Hint{1} + obj.Hint{1}';
%             obj.Hint{1} = obj.Hint{1} + g(1)*kron(obj.JxB,obj.JxS) + g(2)*kron(obj.JyB,obj.JyS) ...
%                 + g(3)*kron(obj.JzB,obj.JzS);
%         end
%         
        function prepareHint(obj,g) 
            %must be called before prepareSim of parent class
            %Hint must be represented in eigenbasis of HS!!!
            %g is a matrix for the interactions: \sum_jk g_jk B_j S_k
            obj.Hint{1} = zeros(obj.dB*obj.dS);
            for j = 1:5
                for k = 1:5
                    obj.Hint{1} = obj.Hint{1} + g(j,k)*kron(obj.JB{j},obj.JS{k});
                end
            end
%             obj.Hint{1} = g(4)*kron(obj.JpB,obj.JpS');
%             obj.Hint{1} = obj.Hint{1} + obj.Hint{1}';
%             obj.Hint{1} = obj.Hint{1} + g(1)*kron(obj.JxB,obj.JxS) + g(2)*kron(obj.JyB,obj.JyS) ...
%                 + g(3)*kron(obj.JzB,obj.JzS);
        end


        function [rhoSt, Jzt, Wt, NjumpsAv, t] = runTheBitch(obj,Nruns,dt,Nt)
            %not parallelized, single core run
            tic
            obj.prepareSim();
            obj.setdtSim(dt);
            %[rhoSt, Jzt, Wt, NjumpsAv, t] = obj.manyRuns(Nruns,Nt,{obj.JzS});
            [rhoSt, ~, Wt, NjumpsAv, t] = obj.manyRuns(Nruns,Nt,{});
            Jzt = zeros(1,1,Nt+1);
            for n=1:obj.dS
                Jzt = Jzt + rhoSt(n,n,:)*obj.JS{3}(n,n);
            end
            Jzt = Jzt(:).';
            toc
        end
        
        function [rhoSt, obst, Wt, NjumpsAv, t] = runTheWorkers(obj,Nruns,dt,Nt,obs)
            %parallelized run
            tic
            obj.prepareSim();
            obj.setdtSim(dt);
            %[rhoSt, Jzt, Wt, NjumpsAv, t] = obj.manyParRuns(Nruns,Nt,{obj.JzS});
            [rhoSt, ~, Wt, NjumpsAv, t] = obj.manyParRuns(Nruns,Nt,{});
            obst = zeros(length(obs),Nt+1);
            
            for n = 1:length(obs)
                for j=1:length(t)
                    obst(n,j) = trace(rhoSt(:,:,j)*obs{n});                
                end
                
            end
            
            toc
        end
        
        function [rhoSt, Jzt, Wt, St, Et, NjumpsAv, t] = runTheBitchNoise(obj,Nruns,dt,Nt,sigma)
            %not parallelized, single core run
            tic
            obj.prepareSim();
            obj.setdtSim(dt);
            %[rhoSt, Jzt, Wt, NjumpsAv, t] = obj.manyRuns(Nruns,Nt,{obj.JzS});
            [rhoSt, ~, Wt, St, Et, NjumpsAv, t] = obj.manyRunsNoise(Nruns,Nt,{},sigma);
            Jzt = zeros(1,1,Nt+1);
            for n=1:obj.dS
                Jzt = Jzt + rhoSt(n,n,:)*obj.JS{3}(n,n);
            end
            Jzt = Jzt(:).';
            toc
        end
        
        function [rhoSt, Jzt, Wt, St, Et, NjumpsAv, t] = runTheWorkersNoise(obj,Nruns,dt,Nt,sigma)
            %parallelized run
            tic
            obj.prepareSim();
            obj.setdtSim(dt);
            %[rhoSt, Jzt, Wt, NjumpsAv, t] = obj.manyParRuns(Nruns,Nt,{obj.JzS});
            [rhoSt, ~, Wt, St, Et, NjumpsAv, t] = obj.manyParRunsNoise(Nruns,Nt,{},sigma);
            Jzt = zeros(1,1,Nt+1);
            for n=1:obj.dS
                Jzt = Jzt + rhoSt(n,n,:)*obj.JS{3}(n,n);
            end
            Jzt = Jzt(:).';
            toc
        end
        
        function getLpm(obj)
        %get Liouvillian based on analytical results for pm interaction of
        %one spin
            dtInt = obj.dtIntB(1);
            g = obj.gsb(4)*2;
            delta = obj.omegaS - obj.omegaB(1);
            lambda = sqrt(delta^2 + g^2);
            A = delta + lambda;
            N = g^2 + A^2;
            k = 2*g*A/N*sin(lambda*dtInt/2);
            C = exp(-1i*(delta+lambda)*dtInt/2)*g^2/N + exp(-1i*(delta-lambda)*dtInt/2)*A^2/N;
            Imc = imag(C);    
            sz = diag([-0.5 0.5]);
            sp = [0 0; 1 0];
            sm = [0 1; 0 0];
            obj.Lpm = 1i*(obj.omegaS+obj.gammaB(1)*Imc)*( obj.rightMultiply(sz) - obj.leftMultiply(sz)) +...
                obj.gammaB(1)*(abs(1-C)).^2* (obj.rightMultiply(sz)*obj.leftMultiply(sz) - 0.25*eye(obj.dS^2)) +...
                obj.gammaB(1)*k^2*obj.pthB{1}(1)*(obj.rightMultiply(sp)*obj.leftMultiply(sm) - 0.5*obj.leftMultiply(sp*sm) - 0.5*obj.rightMultiply(sp*sm)) +...
                obj.gammaB(1)*k^2*obj.pthB{1}(2)*(obj.rightMultiply(sm)*obj.leftMultiply(sp) - 0.5*obj.leftMultiply(sm*sp) - 0.5*obj.rightMultiply(sm*sp));
%             obj.Lpm = 1i*( obj.rightMultiply(sz*(obj.omegaS)) - obj.leftMultiply(sz*(obj.omegaS))) +...
%               obj.gammaB(1)*k^2*obj.pthB{1}(1)*( obj.rightMultiply(sp)*obj.leftMultiply(sm) - ...
%                 0.5*( obj.leftMultiply(sp*sm) + obj.rightMultiply(sp*sm)) ) +...
%                 obj.gammaB(1)*k^2*obj.pthB{1}(2)*( obj.rightMultiply(sm)*obj.leftMultiply(sp) - ...
%                 0.5*( obj.leftMultiply(sm*sp) + obj.rightMultiply(sm*sp)) );
 
        end
        function getLshort(obj,G,gamma)
            %get Liouvillian assuming short time limit, for single spin, no
            %free evolution. Hint = \sum G_k A_k b_k  = \sum A_k B_k
            obj.Lshort = 1i*( obj.rightMultiply(diag(obj.ES)) - obj.leftMultiply(diag(obj.ES)));
            for k = 1:5
                for l = 1:5
                    K(k,l) = trace(diag(obj.pthB{1})*obj.JB{l}'*obj.JB{k})*G(k)*G(l);
                end
            end
            [evec, ev] = eig(K,'vector');
            for j = 1:5
                obj.LOpsShort{end+1} = zeros(obj.dS);
                for k = 1:5
                    obj.LOpsShort{end} = obj.LOpsShort{end} + evec(k,j)*obj.JS{k};
                end
                obj.LOpsShort{end} = sqrt(gamma*ev(j))*obj.LOpsShort{end};
            end
            for j = 1:5
                LL = obj.LOpsShort{j};
                tmp = LL'*LL;
                obj.Lshort = obj.Lshort + obj.rightMultiply(LL')*obj.leftMultiply(LL) - ...
                0.5*( obj.leftMultiply(tmp) + obj.rightMultiply(tmp));
            end 
        end
    end
end
