clear;clc;
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
gammat = linspace(0,1,101);
nbar = linspace(0,10,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;
for i= 1:length(gammat)
     fprintf('(%d)\n',i);
    for j = 1:length(nbar)
        F0(i,j) = findFIopt(0,gammat(i),nbar(j));
        [x1, f1] = myMinSearch(@(theta) findFIopt(theta,gammat(i),nbar(j)),0,pi/2,1e-12);
        [x2, f2] = myMinSearch(@(theta) findFIopt(theta,gammat(i),nbar(j)),pi/2,pi,1e-12);
        if (f2<f1)
            theta(i,j) = x2;
            FI(i,j) = -f2;
        else
            theta(i,j) = x1;
            FI(i,j) = -f1;
        end
    end
end
save('b1F1.mat');