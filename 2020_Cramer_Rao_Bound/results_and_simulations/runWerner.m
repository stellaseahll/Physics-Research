%run four bell states
gt = linspace(0,2*pi,100);
gammat = logspace(-4,1,100);
nbar = 1;
dnbar = 1e-5;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
pw = 0.1:0.1:1;
F1 = cell(length(pw),1);
F2 = cell(length(pw),1);
F3 = cell(length(pw),1);
F4 = cell(length(pw),1);
for k = 1:length(pw)
    for i = 1:length(gt)
        for j = 1:length(gammat)
            s = spSinglePassXYNEW(gt(i),gammat(j),nbar,dnbar,prepareWerner(pw(k)),4,1);
            F = s.getAllFish();
            F1{k}(i,j) = F(1);
            F2{k}(i,j) = F(2);
            F3{k}(i,j) = F(3);
            F4{k}(i,j) = F(4);
        end
    end
    clear s;
    save(sprintf('Werner_nbar%.2f.mat',nbar));
end
