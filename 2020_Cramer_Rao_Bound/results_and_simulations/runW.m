%run W
clear;
gt = linspace(0,2*pi,100);
gammat = logspace(-4,1,100);
nbar = 10;
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = kron([1;0],kron([1;0],[0;1])) + kron([1;0],kron([0;1],[1;0])) +kron([0;1],kron([1;0],[1;0])) ;
rho = rho*rho'/3;

for i = 1:length(gt)
    for j = 1:length(gammat)
        fprintf('(%d,%d)\n',i,j);
        s = spSinglePassXYNEW(gt(i),gammat(j),nbar,dnbar,rho,3,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
        F3(i,j) = F(3);
    end
end
clear s;
save(sprintf('W_nbar%.2f.mat',nbar));

