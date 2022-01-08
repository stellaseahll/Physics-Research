clear; 
np = 12;
rho{1} = [1 0; 0 0];
rho{2} = [0 0; 0 1];
rho{3} = [1 1; 1 1]/2;
filename = sprintf('%dAncillaT2_g_e_p_tol6_gammat0.4.mat',np);
%filename{2} = sprintf('%dAncillaT2_excited.mat',np);
%filename{3} = sprintf('%dAncillaT2_plus.mat',np);
T = 2;
nbar = 1/(exp(1/T)-1);
dnbar = 1e-5;
gammat = [0.4];
gt = [1/100]*pi;
%gammat = [0.01 0.1 1];
%gt = [1/100 1/10 1/4 1/2]*pi;

Fs = cell(3,1);

for n = 1:3
    F = cell(length(gt),length(gammat));
    for i = 1:length(gt)
        for j = 1:length(gammat)
            fprintf('(%i, %i)',i,j);
            s = spSinglePassXY(gt(i),gammat(j),nbar,dnbar*nbar,rho{n},np,1);
            s.maxit = 10000;
            s.tol = 1e-6;
            s.alg = 4;
            F{i,j} = s.getAllFish();
        end
    end
    Fs{n} = F;
end

clear s F
save(filename);

% clear;
% N=500;
% T = 2;
% nbar = 1/(exp(1/T)-1);
% 
% gt = [0.002:0.002:0.498 0.499:0.00002:0.501 0.502:0.002:0.998]*pi;
% % gt = linspace(0,1,N+1)*pi;
% % gt(1) = [];
% % gt(end) = [];
% %singt = (1:N)/N; %linear spacing of plot variable sin(gt)
% %gt = asin(singt);
% 
% %lamb = (1:(N-1))/N; %linear spacing of 1-exp(-gammat*(2*nbar+1))
% gammat = linspace(0,5,N+1);
% gammat(1)=[];
% %gammat = -log(1-lamb)/(2*nbar+1);
% 
% dnbar = 1e-5;
% 
% rhoA = 0.5*[1,1;1,1];
% %rhoA = [1,0;0,0];
% for i = 1:length(gt)
%     i
%     for j = 1:length(gammat)
%         s = spSinglePassXY(gt(i),gammat(j),nbar,dnbar*nbar,rhoA,2,1);
%         s.alg = 3;
%         F = s.getAllFish();
%         F1(i,j) = F(1);
%         F2(i,j) = F(2);
%     end
% end
% % FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1./Fth;
% p1th = nbar/(2*nbar+1);
% FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1/(p1th-p1th^2);
% save('2AncillaT2_gt_gammat_plus.mat');
% %save('1AncillaT0.2_N500_singt_lambda.mat');
% 
% %% next plot
% 
% T = 2;
% nbar = 1/(exp(1/T)-1);
% rhoA = [1,0;0,0];
% 
% for i = 1:length(gt)
%     i
%     for j = 1:length(gammat)
%         s = spSinglePassXY(gt(i),gammat(j),nbar,dnbar*nbar,rhoA,2,1);
%         s.alg = 3;
%         F1(i,j) = s.getAllFish();
%         F = s.getAllFish();
%         F1(i,j) = F(1);
%         F2(i,j) = F(2);
%     end
% end
% %FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1./Fth;
% p1th = nbar/(2*nbar+1);
% FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1/(p1th-p1th^2);
% save('2AncillaT2_gt_gammat.mat');
%save('1AncillaT2_N500_singt_lambda.mat');
