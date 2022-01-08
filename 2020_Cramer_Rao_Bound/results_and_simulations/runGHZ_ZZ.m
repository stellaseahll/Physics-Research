clear;
gt = 0;
gzt = linspace(0,pi/2,100);
gammat = linspace(0,1.5,100);
gammat(1) = 1e-6;
nbar = 1/(exp(1/2)-1);
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = zeros(8);
rho([1 8],[1 8]) = [1 1; 1 1]/2;
for i = 1:length(gammat)
    fprintf('(%d)\n',i);
    for j = 1:length(gzt)

        s = spSinglePassXYZ(gt,gzt(j),gammat(i),nbar,dnbar,rho,3,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
        F3(i,j) = F(3);
    end
end
clear s;
save(sprintf('ZZint_GHZ_nbar%.2f.mat',nbar));

