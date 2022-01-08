clear; 
np = 20;
rho = [1 0];
% rho{3} = [1 1; 1 1]/2;
filename = sprintf('Classical%dAncillaT2_g_e.mat',np);
%filename{2} = sprintf('%dAncillaT2_excited.mat',np);
%filename{3} = sprintf('%dAncillaT2_plus.mat',np);
T = 2;
nbar = 1/(exp(1/T)-1);
dnbar = 1e-5;
gammat = [0.01 0.1 1 10];
gt = [1/100 1/10 1/4 1/2]*pi;

Fs = cell(2,1);

parfor n = 1:2
    F = cell(length(gt),length(gammat));
    for i = 1:length(gt)
        for j = 1:length(gammat)
            for k = 1:np
                F{i,j}(k) =getClassicalFish(gt(i),gammat(j),nbar,dnbar,rho(n),k,1);
            end
        end
    end
    Fs{n} = F;
end

clear F
save(filename);

