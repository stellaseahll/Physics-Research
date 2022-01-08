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
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j)
        [x1, f1] = myMinSearch(@(theta) findFI3optBRUTEFORCE(theta,gammat(i),nbar(j)),0,pi/2,1e-8);
        [x2, f2] = myMinSearch(@(theta) findFI3optBRUTEFORCE(theta,gammat(i),nbar(j)),pi/2,pi,1e-8);
        if (f2<f1)
            theta(i,j) = x2;
            FI(i,j) = -f2;
        else
            theta(i,j) = x1;
            FI(i,j) = -f1;
        end
    end
end
save('optAngleF3.mat');