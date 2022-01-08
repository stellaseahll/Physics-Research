nbar = 1;
gt = (0:0.01:1)*pi;
gt([1 end]) = [];
gammat = logspace(-3,1,100);

for i = 1:length(gt)
    i
    for j = 1:length(gammat)
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 0; 0 0],5,1);
        Fg{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[0 0; 0 1],5,1);
        Fe{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 1; 1 1]/2,5,1);
        Fp{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1+1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1-1/sqrt(2)]/2,5,1);
        F45{i,j} = s.getAllFish(); hold on;
    end
end

save('nbar1.mat');

nbar = 10;
gt = (0:0.025:1)*pi;
gt([1 end]) = [];
gammat = logspace(-3,1,40);

for i = 1:length(gt)
    i
    for j = 1:length(gammat)
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 0; 0 0],5,1);
        Fg{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[0 0; 0 1],5,1);
        Fe{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 1; 1 1]/2,5,1);
        Fp{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1+1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1-1/sqrt(2)]/2,5,1);
        F45{i,j} = s.getAllFish(); hold on;
    end
end

save('nbar10.mat');

nbar = 0.1;
gt = (0:0.01:1)*pi;
gt([1 end]) = [];
gammat = logspace(-3,1,100);

for i = 1:length(gt)
    i
    for j = 1:length(gammat)
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 0; 0 0],5,1);
        Fg{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[0 0; 0 1],5,1);
        Fe{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 1; 1 1]/2,5,1);
        Fp{i,j} = s.getAllFish(); hold on;
        s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1+1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1-1/sqrt(2)]/2,5,1);
        F45{i,j} = s.getAllFish(); hold on;
    end
end

save('nbar0.1.mat');
