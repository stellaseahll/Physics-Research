% Class definition for multipass considering a system ALWAYS
% coupled to a bath, subjected to multipasses by another probe
classdef singlePassXY < handle
    
    properties
        dS = 2; %system dimension
        dP = 2; %probe dimension
        dTot = 4;
        NP = 1; %number of probes
        K = 1; %partial swap parameter, K = sin^2(gtint)
        Lambda = 1; %gamma*twait
        nbar = 1; %system temperature
        dnbar = 0.001;
        Uint %Unitary superop during interaction, includes free Hamiltonian
        L %Dissipation superop
        Mtot %Total superop
        MapSA %Superop on system + 1 ancilla
        JxS
        JyS
        JzS
        JpS
        JmS
        JxP
        JyP
        JzP
        JpP
        JmP
        JSmat %x,y,z,p,m,mp,pm
        JPmat %x,y,z,p,m,mp,pm
        rhoP0 %initial probe state
        rhoSP0 %initial state of S+P
        rhoSS %system steady state
        rhot %before int-afterint-before int-after int....
        DampOps
        gsp
        isSteadOp %decides whether system initially in thermal or steady state
    end
    
    methods
        function obj = singlePassXY(K,lambda,nbar,dnbar,rhoP0,np,isSteadOp)
            obj.nbar = [nbar-dnbar,nbar,nbar+dnbar];
            obj.dnbar = dnbar;
            obj.prepareLambda(lambda);
            obj.prepareMapSA()
            obj.K = K;
            obj.prepareSteadySystem(rhoP0);
            obj.NP = np;
            obj.isSteadOp = isSteadOp;
            obj.dTot = obj.dS*obj.dP^np;
            obj.getOperators();
            obj.prepareProbeState(rhoP0);
            obj.initializeState();
            obj.prepareUnitary();
            obj.prepareDamping();
            obj.loadMap();
            
        end
        
        function obj = prepareLambda(obj,L)
            L1 = 1-(1-L)^((2*obj.nbar(1)+1)/(2*obj.nbar(2)+1));
            L2 = 1-(1-L)^((2*obj.nbar(3)+1)/(2*obj.nbar(2)+1));
            obj.Lambda = [L1 L L2];
        end
        
        function prepareProbeState(obj,rho)
            if (length(rho)== obj.dP)
                obj.rhoP0 = rho ;
                for i = 2:obj.NP
                    obj.rhoP0 = kron(rho,obj.rhoP0);
                end
            elseif (length(rho)== obj.dP^obj.NP)
                obj.rhoP0 = rho;
            else
                error('Check input state of probe');
            end
        end
        
        function prepareUnitary(obj)
            for k = 1:obj.NP
                tmp = -1i*obj.JSmat{4}*obj.JPmat{5,k}*sqrt(obj.K)-1i*obj.JSmat{5}*obj.JPmat{4,k}*sqrt(obj.K)+...
                    obj.JSmat{6}*obj.JPmat{7,k}*sqrt(1-obj.K)+obj.JSmat{7}*obj.JPmat{6,k}*sqrt(1-obj.K)+...
                    obj.JSmat{6}*obj.JPmat{6,k}+obj.JSmat{7}*obj.JPmat{7,k}; %tmp = partial swap
                obj.Uint{k} = obj.lrMultiply(tmp);
            end
        end
        
        function initializeState(obj)
            if (~obj.isSteadOp)
                for i = 1:3
                    obj.rhoSP0{i} = kron(obj.prepareThermalState(obj.nbar(i)),obj.rhoP0);
                end
            else
                for i = 1:3
                    obj.rhoSP0{i} = kron(obj.rhoSS{i},obj.rhoP0);
                end
            end
        end
        
        
        function prepareDamping(obj)
            for i = 1:3
                pg = (obj.nbar(i)+1)/(2*obj.nbar(i)+1);
                obj.L{i} = pg*obj.lrMultiply(obj.DampOps{1,i})...
                    +pg*obj.lrMultiply(obj.DampOps{2,i})...
                    +(1-pg)*obj.lrMultiply(obj.DampOps{3,i})...
                    +(1-pg)*obj.lrMultiply(obj.DampOps{4,i});
            end
            %             n = 1./(exp(obj.beta)+1);
            %             for i = 1:3
            %                 obj.L{i} = sparse(obj.nbar(i)*obj.getDiss(obj.JSmat{4}) +(1+obj.nbar(i))*obj.getDiss(obj.JSmat{5}));
            %             end
        end
        
        function loadMap(obj)
            for i = 1:3
                obj.Mtot{i} = obj.Uint{1};
                for j = 2:obj.NP
                    obj.Mtot{i} = obj.Uint{j}*obj.L{i}*obj.Mtot{i};
                end
            end
        end
        
        function getOperators(obj)
            JS = (obj.dS-1)/2;
            m = (-JS:JS).';
            obj.JzS = diag(-JS:JS);
            t = sqrt((JS-m).*(JS+m+1));
            obj.JpS = diag(t(1:end-1),-1);
            obj.JmS = obj.JpS';
            obj.JxS = 0.5*(obj.JpS + obj.JmS);
            obj.JyS = 0.5*1i*(obj.JpS - obj.JpS');
            obj.JSmat{1} = kron(obj.JxS,eye(obj.dP^obj.NP));
            obj.JSmat{2} = kron(obj.JyS,eye(obj.dP^obj.NP));
            obj.JSmat{3} = kron(obj.JzS,eye(obj.dP^obj.NP));
            obj.JSmat{4} = kron(obj.JpS,eye(obj.dP^obj.NP));
            obj.JSmat{5} = kron(obj.JmS,eye(obj.dP^obj.NP));
            obj.JSmat{6} = kron(obj.JmS*obj.JpS,eye(obj.dP^obj.NP));
            obj.JSmat{7} = kron(obj.JpS*obj.JmS,eye(obj.dP^obj.NP));
            for i = 1:3
                obj.DampOps{1,i} = kron([1 0; 0 sqrt(1-obj.Lambda(i)) ],eye(obj.dP^obj.NP));
                obj.DampOps{2,i} = kron([0 sqrt(obj.Lambda(i)); 0  0],eye(obj.dP^obj.NP));
                obj.DampOps{3,i} = kron([sqrt(1-obj.Lambda(i)) 0; 0 1],eye(obj.dP^obj.NP));
                obj.DampOps{4,i} = kron([0 0; sqrt(obj.Lambda(i)) 0],eye(obj.dP^obj.NP));
            end
            JP = (obj.dP-1)/2;
            obj.JzP = diag(-JP:JP);
            m = (-JP:JP).';
            t = sqrt((JP-m).*(JP+m+1));
            obj.JpP = diag(t(1:end-1),-1);
            obj.JmP = obj.JpP';
            obj.JxP = 0.5*(obj.JpP + obj.JmP);
            obj.JyP = 0.5*1i*(obj.JpP - obj.JpP');
            obj.JPmat = cell(5,obj.NP);
            JpmP = obj.JpP*obj.JmP;
            JmpP = obj.JmP*obj.JpP;
            for j = 1:obj.NP
                obj.JPmat{1,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JxP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{2,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JyP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{3,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JzP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{4,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JpP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{5,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JmP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{6,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(JmpP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{7,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(JpmP,eye(obj.dP^(obj.NP-j))));
            end
        end
        
        function prepareMapSA(obj)
            obj.NP = 1;
            obj.getOperators();
            obj.prepareUnitary();
            obj.prepareDamping();
            obj.MapSA{1} = obj.L{1}*obj.Uint{1};
            obj.MapSA{2} = obj.L{2}*obj.Uint{1};
            obj.MapSA{3} = obj.L{3}*obj.Uint{1};
        end
        
        function prepareSteadySystem(obj,rhoA)
            rho = reshape(kron([1 0; 0 0],rhoA),16,1);
            obj.prepareMapSA();
            for j = 1:3
                for i = 1:500
                    rhos = obj.ptrace(reshape(obj.MapSA{j}*rho,4,4),2,2);
                    rho = reshape(kron(rhos,rhoA),16,1);
                end
                obj.rhoSS{j} = obj.ptrace(reshape(rho,4,4),2,2);
            end

        end
        
%         function prepareSteadySystem(obj,rhoA)
%             CohA = rhoA(1,2);
%             pA = rhoA(1,1);
%             AbsA2 = abs(rhoA(1,2))^2;
%             pg = (obj.nbar+1)/(2*obj.nbar+1);
%             k = obj.K;
%             L = obj.Lambda;
%             denCoh = k.*(1-sqrt(1-k).*sqrt(1-L)+4.*AbsA2.*sqrt(1-k).*sqrt(1-L)).*(-1+L)+(-1+sqrt(1-k).*sqrt(1-L)).*L;
%             numCoh= -CohA.*sqrt(k).*sqrt(1-L).*(L+k.*(-1+L).*(-1+2.*pA)-2.*L.*pg);
%             denp = 4.*(AbsA2).*sqrt(1-k).*k.*(1-L).^(3/2)+k.*(-1+sqrt(1-k).*sqrt(1-L)).*(-1+L)+L-sqrt(1-k).*sqrt(1-L).*L;
%             nump = 2.*AbsA2.*sqrt(1-k).*k.*(1-L).^(3/2)+(-1+sqrt(1-k).*sqrt(1-L)).*(k.*(-1+L).*pA-L.*pg);
%             pSS = nump./denp;
%             CohSS = numCoh./denCoh;
%             for i = 1:3
%                 obj.rhoSS{i} = [pSS(i) CohSS(i); CohSS(i)' 1-pSS(i)];
%             end
%         end

        function rhot = getFinalState(obj,idx)
            rho0 = obj.rhoSP0{idx}(:);
            rhot = obj.Mtot{idx}*rho0;
        end
        
        
        function rho0 = prepareThermalState(obj,nbar)
            rho0 = diag([nbar+1 nbar]/(2*nbar+1));
        end
        
        function F = getFish(obj)
            %get Fisher information before pass Npass
            drho = (obj.getFinalState(3)-obj.getFinalState(1))/2/obj.dnbar;
            rho = obj.getFinalState(2);
            rhoP = obj.getReducedState(rho);
            drhoP = obj.getReducedState(drho);
            F = obj.catchFish(rhoP,drhoP);
        end
        
        function rhoP = getReducedState(obj,rho)
            %get reduced state of probe
            rho = reshape(rho,obj.dTot,obj.dTot);
            %[~,rhoP] = obj.ptrace(rho,obj.dS,obj.dP^obj.NP);
            %can do it directly if tracing over 1st qubit, probably faster
            %for high number of probes
            %WORKS ONLY FOR SYSTEM QUBIT!
            d2 = obj.dP^obj.NP;
            rhoP = rho(1:d2,1:d2) + rho(d2+1:end,d2+1:end);
        end
        
        function F = catchFish(obj,rho,drho)
            drhot = drho.';
            F = 2*drhot(:).'*inv(kron(rho.',eye(length(rho)))+kron(eye(length(rho)),rho))*drho(:);
        end
        
        
        function L = getDiss(obj,A)
            %give the superoperator for the dissipator A
            AdgA = A'*A;
            L = obj.lrMultiply(A) - 0.5*(obj.leftMultiply(AdgA) + obj.rightMultiply(AdgA));
        end
        
        function L = calcFreeEvolve(obj,H)
            %give the superoperator for the dissipator A
            L = -1i* (obj.leftMultiply(H) - obj.rightMultiply(H));
        end
        
        function alpha = getExponent(obj,Fn)
            %takes in a vector Fn and assumrs F_n = k N ^alpha
            N = 1:length(Fn);
            c = polyfit(log(N),log(Fn),1);
            alpha = c(1);
        end
        
        
        function M = leftMultiply(obj,A)
            M = (kron(eye(length(A)),A));
        end
        
        function M = rightMultiply(obj,A)
            M = (kron(A.',eye(length(A))));
        end
        
        function M = lrMultiply(obj,A)
            M = (kron(conj(A),A));
        end
        
        function [A B] = ptrace(obj, AB, dimA, dimB)
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
        
        function I = getInsult(obj)
            insults = {'You suck hard!', 'Your mother stinks like a pile of Wagyu cow shit!', 'You are a worse problem to mankind than climate change!', 'Seeing your face gives me eye cancer.', 'Ever considered a face surgery?', 'Your mother should have considered abortion.', 'Could you just try for a moment to suck a little less, please?', 'Is it your face or your arse I am looking at?','Why dont you play in hell?','You smell from your mouth like a pig from its arse!','Your mother would be proud of you---if she understood anything!','Good! Keep up the great work until you die!','Please turn around. I feel disgusted.','DIVISION ERROR! Tried to divide by your IQ.','OUT OF BOUNDS ERROR in variable yourUgliness!','You again? I think I must puke.','Your fingers just touching my keyboard, what an utterly disgusting sensation!','Dont you think it is time for your monthly wash?','What about taking some holidays? You should visit Shitaly!','You remind me of something. I really need to empty the garbage bin!','The endless stupidity of mankind, how could I forget...','You know what you are doing is completely pointless, right?','C O D T S O P  R U O Y  Y E B O',...
                '________________$$$$\n______________$$____$$\n______________$$____$$\n______________$$____$$\n______________$$____$$\n______________$$____$$\n__________$$$$$$____$$$$$$\n________$$____$$____$$____$$$$\n________$$____$$____$$____$$__$$\n$$$$$$__$$____$$____$$____$$____$$\n$$____$$$$________________$$____$$\n$$______$$______________________$$\n__$$____$$______________________$$\n___$$$__$$______________________$$\n____$$__________________________$$\n_____$$$________________________$$\n______$$______________________$$$\n_______$$$____________________$$\n________$$____________________$$\n_________$$$________________$$$\n__________$$________________$$\n__________$$$$$$$$$$$$$$$$$$$$\n\n',...
                '_______________________##############\n___________________##################\n________________####################\n______________################__####\n____________########___________####\n__________#######___###________####\n________######____#######______###\n______######_____##____####____###\n_____#####______##________##__###\n____####_______##_##___##_##__###\n___####________##__##_##__#######\n__###__________##_________#######\n__##___________###_______########\n__##____________#################\n_______________##################\n_____________################_###\n_____________###############_###\n_____________##############_####\n_____________#############_######\n_____________############_#######\n____________############_########\n____________###########_#########\n___________###########_###########\n___________##########_############\n___________#########_#############\n___________########_##############\n___________#######_###############\n____________#####_###############\n___________#####_################\n__________####_##################\n_________#####_##################\n________#####_###################\n_______####__####################\n______####___####################\n_____####___#####################\n_____####___######################\n_____###____######################\n____________######################\n___________########################\n___________########################\n___________####___########___######\n\n'};
            n = length(insults);
            r = round(1 + rand()*(n-1));
            I = insults{r};
            fprintf(['\n\n',I,'\n\n\n']);
        end
        
        function M = screwYou(obj)
            x = [0.5 0.5 1 1 1 1.5 1.5 1.5 2 2 2 2.5 2.5 2.5 3 3];
            y = [0 1.5 2 1 2.5 2.5 1 4 4 1 2.5 2.5 1 2 2 0];
            X = [];
            Y = [];
            for i = 1:length(x)-1
                X = [X linspace(x(i),x(i+1),10)];
                Y = [Y linspace(y(i),y(i+1),10)];
            end
            for i = 1:length(X)
                plot(X(1:i),Y(1:i),'k');
                axis off;
                axis equal;
                M(i) = getframe;
            end
            for i = length(X)+1:length(X)*2
                if mod(i,5)<=2
                    clf;
                    axis off;
                    axis equal;
                    M(i) = getframe;
                else
                    plot(X,Y,'k');
                    axis off;
                    axis equal;
                    M(i) = getframe;
                end
            end
            
        end
        
    end
    
end


