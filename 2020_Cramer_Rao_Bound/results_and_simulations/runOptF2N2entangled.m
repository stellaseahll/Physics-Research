clf;clear;clc;
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
gammat = linspace(0,1,101);
nbar = linspace(0,10,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;
phi = 0;
options = optimset('MaxFunEvals',1e6,'TolFun',1e-7);
for i= 61:length(gammat)
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        [x1, f1] = fminsearch(@(x) findFI2N2entangled(x,gammat(i),nbar(j)),[pi/8; pi/8;0],options);
        theta1(i,j) = x1(1);
        theta2(i,j) = x1(2);
        phi(i,j) = x1(3);
        F2(i,j) = -f1;
    end
	save('optAngleF2N2entangled_linspace3.mat');
end
