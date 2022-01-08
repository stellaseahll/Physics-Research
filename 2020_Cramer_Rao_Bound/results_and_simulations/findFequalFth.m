clear;clc;
load('F2N2.mat');
nbar_left =  0.01;
nbar_right = 1;
gammat = 1;
tol = 1e-5;
rho = zeros(4,4);
rho(1) = 1;

while (abs(nbar_left-nbar_right)>tol)
    nbar_mid = 0.5*(nbar_left+nbar_right);
    nbar = [nbar_left,nbar_mid,nbar_right];
   
    xguess(1,1) = 0;
    xguess(2,1) = 0;
    xguess(3,1) = 0;
    xguess(4,1) = 0;
    xguess(5,1) = 0;
    options = optimset('MaxFunEvals',1e6,'TolFun',1e-5);
    [x1, fleft] = fminsearch(@(x) findFI2N2var5(x,gammat,nbar_left),xguess,options);
    [x2, fright] = fminsearch(@(x) findFI2N2var5(x,gammat,nbar_right),xguess,options);
    [x3, fmid] = fminsearch(@(x) findFI2N2var5(x,gammat,nbar_mid),xguess,options);
    Fth = 2./(1+2*nbar).^2./(1+nbar)./nbar
    FI = [-fleft,-fmid,-fright]
    Fdiff = sign(FI-Fth)
    if (Fdiff(1)==Fdiff(2))
        nbar_left = nbar_mid;
    else
        nbar_right = nbar_mid;
    end
    nbar_left
    nbar_right
end