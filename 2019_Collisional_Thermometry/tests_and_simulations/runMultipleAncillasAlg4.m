clear;
gt = 1/100*pi;
for N = 14%:15
    for state = 1:3
        load(sprintf('%dAncillaT2_State%d_gt%.2f.mat',N-1,state,gt/pi));
        filename = sprintf('%dAncillaT2_State%d_gt%.2f.mat',N,state,gt/pi);
        N=14;
        for j = 1:length(gammat)
            fprintf('(%d,%d,%d)',N,state,j);
            s = spSinglePassXYNew(gt,gammat(j),nbar,dnbar*nbar,rho{state},N,1);
            s.alg = 4;
            F(j,N) = s.getQuickFish();
            clear s;
            save(filename);
        end
    end
end