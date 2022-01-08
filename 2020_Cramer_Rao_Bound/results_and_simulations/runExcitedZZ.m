%run ground state
clear;
gt = 0;
gzt = linspace(0,pi,100);
gammat = logspace(-4,0.3,100);
gammat(end) = 2;
nbar = 1;
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = [1 1; 1 1]/2;
for i = 1:length(gammat)
    fprintf('(%d)\n',i);
    for j = 1:length(gzt)
        s = spSinglePassXYZ(gt,gzt(j),gammat(i),nbar,dnbar,rho,4,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
        F3(i,j) = F(3);
        F4(i,j) = F(4);
    end
end
clear s;
save(sprintf('ZZint_plus%.2f.mat',nbar));

