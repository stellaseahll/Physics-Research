%run GHZ
clear;
gt = linspace(0,2*pi,100);
gammat = logspace(-4,1,100);
nbar = 1;
dnbar = 1e-5;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = zeros(8);
rho([1 8],[1 8]) = [1 1; 1 1]/2;

for i = 1:length(gt)
    for j = 1:length(gammat)
        s = spSinglePassXYNEW(gt(i),gammat(j),nbar,dnbar,rho,3,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
        F3(i,j) = F(3);
    end
end
clear s;
save(sprintf('GHZ_nbar%.2f.mat',nbar));

