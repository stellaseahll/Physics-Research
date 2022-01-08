clear;clc;
load('F2N2.mat','uncorr_F2','uncorr_phi','uncorr_theta1','uncorr_theta2');
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
gammat = linspace(0,1,101);
nbar = logspace(-3,0,101);
gammat(1) = 0.0001;
% nbar(1) = 0.001;
options = optimset('MaxFunEvals',1e6,'TolFun',1e-6);
xinput = [0 0 0 0 0; 0 pi/2 0 0 0; pi/2 0 0 0 0]';
for i= 1:length(gammat)
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
%         options = optimset('MaxFunEvals',1e6,'TolFun',10^min([floor(log10(uncorr_F2(i,j)))-4,-6]));
        fgg = findFI2N2var5([0 0 0 0 0],gammat(i),nbar(j));
        fgp = findFI2N2var5([0 pi/2 0 0 0],gammat(i),nbar(j));
        fpg = findFI2N2var5([pi/2 0 0 0 0],gammat(i),nbar(j));
        f = -[fgg fgp fpg];
        idx = find(f==max(f));
        f = max(f);
        options = optimset('MaxFunEvals',1e6,'TolFun',10^min([floor(log10(f))-4,-6]));
        [x1, f1] = fminsearch(@(x) findFI2N2var5(x(:,idx),gammat(i),nbar(j)),xinput(:,idx),options);
        while (-f1<f && it<5)
            it=it+1;
            xguess = xinput(:,idx) + (rand(5,1)-0.5)*0.2;
            [x1, f1] = fminsearch(@(x) findFI2N2var5(x,gammat(i),nbar(j)),xguess,options);   
        end
        if (-f1>f)
            b2_theta1(i,j) = x1(1);
            b2_theta2(i,j) = x1(2);
            b2_phi(i,j) = x1(3);
            b2_alpha(i,j) = x1(4);
            b2_coeff(i,j) = x1(5);
            b2_F2(i,j) = -f1;
        else
            b2_theta1(i,j) = xinput(1,idx);
            b2_theta2(i,j) = xinput(2,idx)
            b2_phi(i,j) = xinput(3,idx)
            b2_alpha(i,j) = xinput(4,idx)
            b2_coeff(i,j) = xinput(5,idx)
            b2_F2(i,j) = f;
        end
        if (b2_F2(i,j)-f<0)
            fprintf('error');
        end
    end
  save('optF2N2_logspace.mat');
end