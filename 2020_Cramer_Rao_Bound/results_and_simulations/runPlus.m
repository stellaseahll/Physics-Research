function [gt gammat F1 F2 Fth]=runPlus(nbar)

gt = linspace(0,pi,102);
gt(1) = [];
gt(end) = [];
gammat = logspace(-2,1,100);
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = [1 0; 0 0];
for i = 1:length(gt)
    fprintf('(%d)\n',i);
    for j = 1:length(gammat)
        
        s = spSinglePassXYNEW(gt(i),gammat(j),nbar,dnbar,rho,2,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
    end
end
clear s;
