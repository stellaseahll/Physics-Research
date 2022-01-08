clear;
gt = pi/100;
gammat = 0.1;
nbar = logspace(-2,1,100);
dnbar = 1e-5;
np = 7;

for i = 1:length(nbar)
    sg = spSinglePassXY(gt,gammat,nbar(i),dnbar*nbar(i),[1 0; 0 0],np,1);
    Fg(i,:) = sg.getAllFish();
    se = spSinglePassXY(gt,gammat,nbar(i),dnbar*nbar(i),[0 0; 0 1],np,1);
    Fe(i,:) = se.getAllFish();
end

save('data2.mat');