% Parameters
dimS = 2;
omegaS = 10;
dimB = 2;
omegaB = 10;
TB = 10;
gamma = 1e-2;
gxx = 1;
gyy = 1;
gzz = 1; 
gpm = 0;
gsb = [gxx gyy gzz gpm];
Ntrials = 10000;
dtInt = 1;
sigma = [0.01 0.1 1];
rhoSt = cell(length(dtInt),length(sigma));
rhoTh = exp(-(1:dimS)*omegaS/TB);
rhoTh = diag(rhoTh)/sum(rhoTh);

for m = 1:length(dtInt)
    for n = 1:length(sigma)
        rhoSt{m,n} = zeros(dimS);
        wB = (1+(rand(1,Ntrials)*2*sigma(n)-sigma(n)))*omegaB;
        for j = 1:Ntrials
            if ~mod(j,1000)
               fprintf('(%d,%d,%d)\n',m,n,j);
            end            
            c = collisionModelSpinBath(dimS,omegaS,dimB,wB(j),TB,gamma,dtInt(m),gsb,[1 0; 0 0]);
            c.findSSwithoutSim();
            rhoSt{m,n} = rhoSt{m,n}+c.rhoSSscatter{1}/Ntrials;
        end
    end
end
% dtInt = logspace(-2,2,200);
% wB = [0.1 10 20];
% for j = 1:length(wB)
%     j
%     for k = 1:length(dtInt)
%         c = collisionModelSpinBath(dimS,omegaS,dimB,wB(j),TB,gamma,dtInt(k),gsb,[1 0; 0 0]);
%         c.findSSwithoutSim();
%         rhoSt{j,k} = c.rhoSSscatter{1};
%         p0(j,k) = rhoSt{j,k}(1);
%         pthB(j,k) = c.pthB{1}(1);
%     end
% end
%     
%     
