function [F1 F2]=runPlusXYZ(nbar,gammat)

gt = linspace(0,2*pi,100);
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = [1 1; 1 1]/2;
for i = 1:length(gt)
    fprintf('(%d)\n',i);
    for j = 1:length(gzt)
        
        s = spSinglePassXYZ(gt(i),gzt(j),gammat,nbar,dnbar,rho,2,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
    end
end
clear s;
save(sprintf('XYZint_plus_nbar%.2f_gammat_%.3f.mat',nbar,gammat));

