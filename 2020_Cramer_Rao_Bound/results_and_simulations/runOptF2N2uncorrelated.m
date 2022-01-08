clf;clear;clc;
load('optF2N2.mat');
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
gammat = linspace(0,1,101);
nbar = linspace(0,10,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;
options = optimset('MaxFunEvals',1e6,'TolFun',1e-6);
for i= 76:101
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        options = optimset('MaxFunEvals',1e6,'TolFun',10^min([floor(log10(F2N2b2(i,j)))-4,-5]));
        switch region(i,j)
            case 1
                [x1, f1] = fminsearch(@(x) findFI2N2uncorr(x,gammat(i),nbar(j)),[0.1; pi/2.1;0.01],options);
            case 2
                [x1, f1] = fminsearch(@(x) findFI2N2uncorr(x,gammat(i),nbar(j)),[0.1; 0.1;0.01],options);
            case 3
                [x1, f1] = fminsearch(@(x) findFI2N2uncorr(x,gammat(i),nbar(j)),[pi/2;0.1;0.01],options);
        end
        theta1uncorr(i,j) = x1(1);
        theta2uncorr(i,j) = x1(2);
        phiuncorr(i,j) = x1(3);
        F2uncorr(i,j) = max([-f1 F2N2b2(i,j)]);
    end
	save('optAngleF2N2uncorrelated_linspace3.mat');
end
