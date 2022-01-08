%run four bell states
gt = linspace(0,2*pi,100);
gammat = logspace(-4,1,100);
nbar = 0.1;
dnbar = 1e-5;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
F1 = cell(4,1);
F2 = cell(4,1);
F3 = cell(4,1);
F4 = cell(4,1);
for k = 1:4
    for i = 1:length(gt)
        for j = 1:length(gammat)
            fprintf('(%d,%d)\n',i,j);
            s = spSinglePassXYNEW(gt(i),gammat(j),nbar,dnbar,prepareBell(k),4,1);
            F = s.getAllFish();
            F1{k}(i,j) = F(1);
            F2{k}(i,j) = F(2);
            F3{k}(i,j) = F(3);
            F4{k}(i,j) = F(4);
        end
    end
    clear s;
    save(sprintf('Bell_nbar%.2f.mat',nbar));
end
