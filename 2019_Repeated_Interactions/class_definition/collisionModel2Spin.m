%Class definition for repeated interaction at Poisson times, system and
%multiple baths. System is 2 spins, 2 baths

classdef collisionModel2Spin < collisionModel
    
    properties 
        %additional to parent class
        %spin operators for system 1
        JzS1
        JxS1
        JyS1
        JpS1
        %spin operators for system 2
        JzS2
        JxS2
        JyS2
        JpS2
        %JmS
        %spin operators for bath 1
        JzB1
        JxB1
        JyB1
        JpB1
        %JmB
        %spin operators for bath 2
        JzB2
        JxB2
        JyB2
        JpB2
        %coupling coefficients
        %gxx = 0
        %gyy = 0;
        %gzz = 0;
        %gpm = 0;
    end
    
    methods
        function obj = collisionModel2Spin(ds1,omegas1,ds2,omegas2,db,omegab,Hb1,Hb2,Tb,gamma,dtint,gsb1,gsb2,gss,rho0)
            % g matrix = [gxx gyy gzz gpm]
            % Hb = Hamiltonian of bath
            % rho0 in computational basis
            %call constructor of parent class
            obj = obj@collisionModel(ds1*ds2,2,db,omegab,Tb,gamma,dtint);
            IdS1 = eye(ds1);
            IdS2 = eye(ds2);
            %operators for system 1 acting on S1 tensor S2
            J = (ds1-1)/2;
            m = (-J:J).';
            obj.JzS1 = kron(diag(m),IdS2);
            t = sqrt((J-m).*(J+m+1));
            obj.JpS1 = kron(diag(t(1:end-1),-1),IdS2);
            obj.JxS1 = 0.5*(obj.JpS1 + obj.JpS1');
            obj.JyS1 = 0.5/1i*(obj.JpS1 - obj.JpS1'); %CHECK!
            
            %operators for system 2 acting on S1 tensor S2
            J = (ds2-1)/2;
            m = (-J:J).';
            obj.JzS2 = kron(IdS1,diag(m));
            t = sqrt((J-m).*(J+m+1));
            obj.JpS2 = kron(IdS1,diag(t(1:end-1),-1));
            obj.JxS2 = 0.5*(obj.JpS2 + obj.JpS2');
            obj.JyS2 = 0.5/1i*(obj.JpS2 - obj.JpS2'); %CHECK!
            
            %prepare system operators, system in ground state (computation)
            obj.prepareHS(gss,omegas1,omegas2)
            obj.rhoS0 = rho0;
            obj.prepareSystemOps()
            %spin operators S

            %spin operators B1
            J = (db(1)-1)/2;
            m = (-J:J).';
            obj.JzB1 = diag(m);
            t = sqrt((J-m).*(J+m+1));
            obj.JpB1 = diag(t(1:end-1),-1);
            obj.JxB1 = 0.5*(obj.JpB1 + obj.JpB1');
            obj.JyB1 = 0.5*1i*(obj.JpB1 - obj.JpB1');
            obj.prepareHb1(Hb1);
            obj.rhoB{1} = expm(-obj.HB{1}/obj.TB(1));
            obj.rhoB{1} = obj.rhoB{1}/trace(obj.rhoB{1});
            %spin operators B2
            J = (db(2)-1)/2;
            m = (-J:J).';
            obj.JzB2 = diag(m);
            t = sqrt((J-m).*(J+m+1));
            obj.JpB2 = diag(t(1:end-1),-1);
            obj.JxB2 = 0.5*(obj.JpB2 + obj.JpB2');
            obj.JyB2 = 0.5*1i*(obj.JpB2 - obj.JpB2');
            obj.prepareHb2(Hb2);
            obj.rhoB{2} = expm(-obj.HB{2}/obj.TB(2));
            obj.rhoB{2} = obj.rhoB{2}/trace(obj.rhoB{2});
            %prepare interaction Hamiltonianc
            obj.prepareHint1(gsb1);
            obj.prepareHint2(gsb2);
% 			obj.getInsult();
        end
        
        function prepareHS(obj,g,omegas1,omegas2)
            %must be called before prepareSim of parent class
            %Hint must be represented in eigenbasis of HS!!!
            obj.HS = g(4)*obj.JpS1*obj.JpS2';
            obj.HS = obj.HS + obj.HS';
            obj.HS = obj.HS + g(1)*obj.JxS1*obj.JxS2 + g(2)*obj.JyS1*obj.JyS2  ...
                + g(3)*obj.JzS1*obj.JzS2 ;
            obj.HS = obj.HS + omegas1*obj.JzS1 + omegas2*obj.JzS2;
        end
        
        function prepareSystemOps(obj)
            [V, h] = eig(obj.HS,'vector'); %A*V = V*D
            obj.Ucomp2E = V;
            obj.ES = h;
            obj.rhoS0 = obj.basisTransform(obj.rhoS0);
        end
        
        function X = basisTransform(obj,X)
            X = obj.Ucomp2E*X*obj.Ucomp2E';
        end
        
        function prepareHb1(obj,Hb) 
            %must be called before prepareSim of parent class
            %Hint must be represented in eigenbasis of HS!!!
            obj.HB{1} = Hb(1)*obj.JxB1 + Hb(2)*obj.JyB1 + Hb(3)*obj.JzB1;
        end
        
        function prepareHb2(obj,Hb) 
            %must be called before prepareSim of parent class
            %Hint must be represented in eigenbasis of HS!!!
            obj.HB{2} = Hb(1)*obj.JxB2 + Hb(2)*obj.JyB2 + Hb(3)*obj.JzB2;
        end
        
        function prepareHint1(obj,g) 
            %must be called before prepareSim of parent class
            %Hint must be represented in eigenbasis of HS!!!
            jxS1 = obj.basisTransform(obj.JxS1);
            jyS1 = obj.basisTransform(obj.JyS1);
            jzS1 = obj.basisTransform(obj.JzS1);
            jpS1 = obj.basisTransform(obj.JpS1);
            
            obj.Hint{1} = g(4)*kron(obj.JpB1,jpS1');
            obj.Hint{1} = obj.Hint{1} + obj.Hint{1}';
            obj.Hint{1} = obj.Hint{1} + g(1)*kron(obj.JxB1,jxS1) + g(2)*kron(obj.JyB1,jyS1) ...
                + g(3)*kron(obj.JzB1,jzS1);
        end
        
        function prepareHint2(obj,g) 
            %must be called before prepareSim of parent class
            %Hint must be represented in eigenbasis of HS!!!
            jxS2 = obj.basisTransform(obj.JxS2);
            jyS2 = obj.basisTransform(obj.JyS2);
            jzS2 = obj.basisTransform(obj.JzS2);
            jpS2 = obj.basisTransform(obj.JpS2);
            
            obj.Hint{2} = g(4)*kron(obj.JpB2,jpS2');
            obj.Hint{2} = obj.Hint{2} + obj.Hint{2}';
            obj.Hint{2} = obj.Hint{2} + g(1)*kron(obj.JxB2,jxS2) + g(2)*kron(obj.JyB2,jyS2) ...
                + g(3)*kron(obj.JzB2,jzS2);
        end

        function [rhoSt, Jzt, NjumpsAv, t] = runTheBitch(obj,Nruns,dt,Nt)
            %not parallelized, single core run
            tic
            obj.prepareSim();
            obj.setdtSim(dt);
            obsv{1} = obj.basisTransform(obj.JzS1);
            obsv{2} = obj.basisTransform(obj.JzS2);
            obsv{3} = diag(obj.ES);
            [rhoSt, Jzt, NjumpsAv, t] = obj.manyRuns(Nruns,Nt,obsv);
            toc
        end
        
        function [rhoSt, obst, NjumpsAv, t] = runTheWorkers(obj,Nruns,dt,Nt)
            %parallelized run
            %rhoSt is in energy eigenbasis.
            tic
            obj.prepareSim();
            obj.setdtSim(dt);
            obs{1} = obj.basisTransform(obj.JzS1);
            obs{2} = obj.basisTransform(obj.JzS2);
            obs{3} = obj.basisTransform(obj.JxS1);
            obs{4} = obj.basisTransform(obj.JxS2);
            obs{5} = diag(obj.ES);
            [rhoSt, ~, Wt, NjumpsAv, t] = obj.manyParRuns(Nruns,Nt,obs);
            for n = 1:length(obs)
                for j=1:length(t)
                    obst(n,j) = trace(rhoSt(:,:,j)*obs{n});                
                end
                
            end
            toc
        end
        
    end
end
