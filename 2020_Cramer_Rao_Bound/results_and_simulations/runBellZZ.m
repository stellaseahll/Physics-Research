%run four bell states
clear;
gt = 0;
gzt = linspace(0,pi/2,100);
gammat = linspace(0,2,100);
gammat(1) = 1e-6;
nbar = 1/(exp(1/2)-1);
dnbar = 1e-5*nbar;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
F1 = cell(4,1);
F2 = cell(4,1);
F3 = cell(4,1);
F4 = cell(4,1);
for k = 1:4
    for i = 1:length(gammat)
        fprintf('(%d,%d)\n',k,i);
        for j = 1:length(gzt)
            s = spSinglePassXYZ(gt,gzt(j),gammat(i),nbar,dnbar,prepareBell(k),4,1);
            F = s.getAllFish();
            F1{k}(i,j) = F(1);
            F2{k}(i,j) = F(2);
            F3{k}(i,j) = F(3);
            F4{k}(i,j) = F(4);
        end
    end
    clear s;
    save(sprintf('ZZint_Bell_nbar%.2f.mat',nbar));
end
