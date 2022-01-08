%run ground state
clear;
gt = 0;
gzt = pi/4;
gammat = logspace(-4,0,100);
nbar = logspace(-3,1,100);
dnbar = 1e-5*nbar;

rho = [1 1; 1 1]/2;
for i = 1:length(gammat)
    fprintf('(%d)\n',i);
    for j = 1:length(nbar)
        Fth(i,j) = 1/(1+2*nbar(j))^2/(1+nbar(j))/nbar(j);
        s = spSinglePassXYZ(gt,gzt,gammat(i),nbar(j),dnbar(j),rho,2,1);
%         s.alg = 4;
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
    end
end
clear s;
% save(sprintf('ZZint_plus_nbar%.2f.mat',nbar));

