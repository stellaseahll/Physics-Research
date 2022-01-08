% Class definition for multipass considering a system ALWAYS
% coupled to a bath, subjected to multipasses by another probe
classdef singlePassMultiSpin < handle
    
    properties
        dS = 2; %system dimension
        dP = 2; %probe dimension
        dTot = 4;
        NP = 1; %number of probes
        omegaP = 1; %array of bath gaps, length NP
        K = 1; %partial swap parameter, K = cos^2(gtint)
        Gammat = 1; %gamma*twait
        omegaS = 1; %system frequency
        TS = 1; %system temperature
        dTS = 0.001;
        Uint %Unitary during interaction, includes free Hamiltonian
        L
        L0
        Lint
        M0  %maps
        Mint
        Mtot
        Mint0 %interaction followed by free
        M0int %free followed by interaction
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
        JSmat %x,y,z,p,m
        JPmat %x,y,z,p,m
        rhoP0 %initial probe state
        rhoSP0 %initial state of S+P
        rhot %before int-afterint-before int-after int....
        gsp
    end
    
    methods
        function obj = singlePassMultiSpin(np,gsp,K,Gammat,rhoP0,Ts,dTs)
            obj.tInt = tInt;
            obj.tWait = tWait;
            obj.omegaS = os;
            obj.omegaP = op;
            obj.NP = np;
            obj.prepareProbeState(rhoP0);
            obj.gamma = gamma;
            obj.dTot = obj.dS*obj.dP^np;
            obj.TS = [Ts-dTs,Ts,Ts+dTs];
            obj.dTS = dTs;
            obj.getOperators();
            obj.loadHamiltonian(gsp);
            obj.loadFreeEvolve();
            obj.loadDissipator();
            obj.loadMap();
            obj.initializeState();
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
        
        function initializeState(obj)
            for i = 1:3
                obj.rhoSP0{i} = kron(obj.prepareThermalState(obj.TS(i)),obj.rhoP0);
            end
        end
        
        function loadHamiltonian(obj,gsp)
            obj.H0 = obj.omegaS*obj.JSmat{3};
            for i = 1:obj.NP
                obj.H0 = obj.H0 + obj.omegaP*obj.JPmat{3,i};
            end
            obj.Hint = cell(1,obj.NP);
            for i = 1:obj.NP
                obj.Hint{i} = zeros(obj.dTot);
                for j = 1:5
                    for k = 1:5
                        obj.Hint{i} = obj.Hint{i} + gsp(j,k)*obj.JSmat{j}*obj.JPmat{k,i};
                    end
                end
            end
        end
        
        function loadFreeEvolve(obj)
            obj.U0 = sparse(obj.calcFreeEvolve(obj.H0));
            obj.Uint = cell(1,obj.NP);
            for i = 1:obj.NP
                obj.Uint{i} = sparse(obj.calcFreeEvolve(obj.Hint{i}));
            end
        end
        
        function loadDissipator(obj)
            n = 1./(exp(obj.omegaS./obj.TS)-1);
            for i = 1:3
                obj.L{i} = sparse(obj.gamma*n(i)*obj.getDiss(obj.JSmat{4}) + obj.gamma*(n(i)+1)*obj.getDiss(obj.JSmat{5}));
                obj.L0{i} = obj.L{i}+obj.U0;
                for j = 1:obj.NP
                    obj.Lint{i,j} = obj.L{i}+obj.Uint{j};
                end
            end
        end
        
        function loadMap(obj)
            for i = 1:3
                obj.M0{i} = sparse(expm(obj.L0{i}*obj.tWait));
                for j = 1:obj.NP
                    obj.Mint{i,j} = sparse(expm(obj.Lint{i,j}*obj.tInt));
                end
            end
            for i = 1:3
                obj.Mtot{i} = obj.Mint{i,1};
                for j = 2:obj.NP
                    obj.Mtot{i} = obj.Mint{i,j}*obj.M0{i}*obj.Mtot{i};
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
            JP = (obj.dP-1)/2;
            obj.JzP = diag(-JP:JP);
            m = (-JP:JP).';
            t = sqrt((JP-m).*(JP+m+1));
            obj.JpP = diag(t(1:end-1),-1);
            obj.JmP = obj.JpP';
            obj.JxP = 0.5*(obj.JpP + obj.JmP);
            obj.JyP = 0.5*1i*(obj.JpP - obj.JpP');
            obj.JPmat = cell(5,obj.NP);
            for j = 1:obj.NP
                obj.JPmat{1,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JxP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{2,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JyP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{3,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JzP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{4,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JpP,eye(obj.dP^(obj.NP-j))));
                obj.JPmat{5,j} = kron(eye(obj.dS*obj.dP^(j-1)),kron(obj.JmP,eye(obj.dP^(obj.NP-j))));
            end
        end
        
        function rhot = getFinalState(obj,idx)
            rho0 = obj.rhoSP0{idx}(:);
            rhot = obj.Mtot{idx}*rho0;
        end
        
        
        function rho0 = prepareThermalState(obj,T)
            pth = exp(-obj.omegaS/T * (0:(obj.dS-1)).' );
            pth = pth/sum(pth);
            rho0 = diag(pth);
        end
        
        function F = getFish(obj)
            %get Fisher information before pass Npass
            drho = (obj.getFinalState(3)-obj.getFinalState(1))/2/obj.dTS;
            rho = obj.getFinalState(2);
            rhoP = obj.getReducedState(rho);
            drhoP = obj.getReducedState(drho);
            F = obj.catchFish(rhoP,drhoP);
        end
        
        function rhoP = getReducedState(obj,rho)
            %get reduced state of probe
            rho = reshape(rho,obj.dTot,obj.dTot);
            [~,rhoP] = obj.ptrace(rho,obj.dS,obj.dP^obj.NP);
        end
        
        function F = catchFish(obj,rho,drho)
            drhot = drho.';
            F = 2*drhot(:).'/inv(kron(rho.',eye(length(rho)))+kron(eye(length(rho)),rho))*drho(:);
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


