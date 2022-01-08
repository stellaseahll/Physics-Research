clf;clear;clc;
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
load('optF2N2.mat');
gammat = linspace(0,1,101);
nbar = linspace(0,10,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;
phi = 0;
options = optimset('MaxFunEvals',1e6,'TolFun',1e-7);
for i= 52%:length(gammat)
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        xguess(1,1) = theta1maxEnt(i,j);
        xguess(2,1) = theta2maxEnt(i,j);
        xguess(3,1) = 0;
        xguess(4,1) = phimaxEnt(i,j);
        options = optimset('MaxFunEvals',1e6,'TolFun',10^min([floor(log10(F2maxEnt(i,j)))-4,-6]));
        [x1, f1] = fminsearch(@(x) findFI2N2entangledvar4(x,gammat(i),nbar(j)),xguess,options);        
        it = 0;
        while (-f1<F2maxEnt(i,j) && it<5)
            it=it+1;
            xguess(1,1) = theta1maxEnt(i,j)+rand*0.2;
            xguess(2,1) = theta2maxEnt(i,j)+rand*0.2;
            xguess(3,1) = rand*0.2;
            xguess(4,1) = phimaxEnt(i,j)+rand*0.2;
           
            [x1, f1] = fminsearch(@(x) findFI2N2entangledvar4(x,gammat(i),nbar(j)),xguess,options);
        end
        if (-f1<F2maxEnt(i,j))
            theta1Ent(i,j) = x1(1);
            theta2Ent(i,j) = x1(2);
            phivarEnt(i,j) = x1(3);
            alphaEnt(i,j) = x1(4);
            F2entvar4(i,j) = -f1;
        else
            theta1Ent(i,j) = theta1maxEnt(i,j);
            theta2Ent(i,j) = theta2maxEnt(i,j);
            phivarEnt(i,j) = 0;
            alphaEnt(i,j) = phimaxEnt(i,j);
            F2entvar4(i,j) = F2maxEnt(i,j);
        end
    end
    save('optAngleF2N2entangledvar4_linspace1.mat');
end
