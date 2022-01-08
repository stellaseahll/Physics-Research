clear;
N=500;
T = 2;
nbar = 1/(exp(1/T)-1);

gt = [0.002:0.002:0.498 0.499:0.00002:0.501 0.502:0.002:0.998]*pi;
% gt = linspace(0,1,N+1)*pi;
% gt(1) = [];
% gt(end) = [];
%singt = (1:N)/N; %linear spacing of plot variable sin(gt)
%gt = asin(singt);

%lamb = (1:(N-1))/N; %linear spacing of 1-exp(-gammat*(2*nbar+1))
gammat = linspace(0,5,N+1);
gammat(1)=[];
%gammat = -log(1-lamb)/(2*nbar+1);

dnbar = 1e-5;

rhoA = 0.5*[1,1;1,1];
%rhoA = [1,0;0,0];
for i = 1:length(gt)
    i
    for j = 1:length(gammat)
        s = spSinglePassXY(gt(i),gammat(j),nbar,dnbar*nbar,rhoA,1,1);
        s.alg = 3;
        F1(i,j) = s.getAllFish();
        %rhoP = s.getReducedState(s.rhot{2});
%        Fth(i,j) = (rhoP(end)- rhoP(end)^2);
%         F1(i,j) = F{i,j}(1);
%         F2(i,j) = F{i,j}(2);
%         F3(i,j) = F{i,j}(3);
%         F4(i,j) = F{i,j}(4);
%         F5(i,j) = F{i,j}(5);
%         F6(i,j) = F{i,j}(6);
%         F7(i,j) = F{i,j}(7);
    end
end
% FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1./Fth;
p1th = nbar/(2*nbar+1);
FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1/(p1th-p1th^2);
save('1AncillaT2_N500_gt_gammat_plus.mat');
%save('1AncillaT0.2_N500_singt_lambda.mat');

%% next plot

T = 2;
nbar = 1/(exp(1/T)-1);
rhoA = [1,0;0,0];

for i = 1:length(gt)
    i
    for j = 1:length(gammat)
        s = spSinglePassXY(gt(i),gammat(j),nbar,dnbar*nbar,rhoA,1,1);
        s.alg = 3;
        F1(i,j) = s.getAllFish();
        %rhoP = s.getReducedState(s.rhot{2});
        %Fth(i,j) = (rhoP(end)- rhoP(end)^2);
%         F1(i,j) = F{i,j}(1);
%         F2(i,j) = F{i,j}(2);
%         F3(i,j) = F{i,j}(3);
%         F4(i,j) = F{i,j}(4);
%         F5(i,j) = F{i,j}(5);
%         F6(i,j) = F{i,j}(6);
%         F7(i,j) = F{i,j}(7);
    end
end
%FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1./Fth;
p1th = nbar/(2*nbar+1);
FFthRat = (exp(1/T)/(exp(1/T)-1)^2)^2*F1/(p1th-p1th^2);
save('1AncillaT2_N500_gt_gammat.mat');
%save('1AncillaT2_N500_singt_lambda.mat');
