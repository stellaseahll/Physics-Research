%run ground state
clear;
gt = linspace(0,2*pi,100);
gzt = linspace(0,2*pi,100);
gammat = 1;
nbar = 1;
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
rho = [1 0; 0 0];
for i = 1:length(gt)
    fprintf('(%d)\n',i);
    for j = 1:length(gzt)
        
        s = spSinglePassXYZ(gt(i),gzt(j),gammat,nbar,dnbar,rho,4,1);
        F = s.getAllFish();
        F1(i,j) = F(1);
        F2(i,j) = F(2);
        F3(i,j) = F(3);
        F4(i,j) = F(4);
    end
end
clear s;
save(sprintf('XYZint_ground_nbar%.2f.mat',nbar));

