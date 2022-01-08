% run ground state
clear;
gt = pi/2;
Gammat = 1;
nbar = logspace(-2,2,50);
dnbar = 1e-6*nbar;
gammat = Gammat./(2*nbar+1);
Fth = 1./(1+2*nbar).^2./(1+nbar)./nbar;
theta = linspace(0,pi/2,50);
for i = 1:length(nbar)
    fprintf('(%d)\n',i);
    for j = 1:length(theta)
        rho = [cos(theta(j)); sin(theta(j))]*[cos(theta(j)); sin(theta(j))]';
        s = spSinglePassXYNEW(gt,gammat(i),nbar(i),dnbar(i),rho,2,1);
        F = s.getAllFish();
        s.getRhoA;
        s.getRhoAA;
        MI(i,j) = s.getMI();
        F1(i,j) = F(1);
        F2(i,j) = max(F(2),F(1)*2);
    end
end
subplot(1,2,1); hold on;
title('Mutual info');
xlabel('theta/pi');
ylabel('log(nbar)');
contourf(theta/pi,log10(nbar),MI);
colorbar;
subplot(1,2,2); hold on;
title('F2/2F1');
xlabel('theta/pi');
ylabel('log(nbar)');
contourf(theta/pi,log10(nbar),real(F2./F1/2));
colorbar;
% clear s;
% % save(sprintf('ground_nbar%.3f.mat',nbar));
% 
% clear;
% gt = pi/2;
% Nt = 1000;
% gammat = rand(1,Nt)*10;
% nbar = rand(1,Nt)*10;
% dnbar = 1e-6*nbar;
% Fth = 1./(1+2*nbar).^2./(1+nbar)./nbar;
% theta = rand(1,Nt)*pi/2;
% for i = 1:Nt
%     fprintf('(%d)\n',i);
%     rho = [cos(theta(i)); sin(theta(i))]*[cos(theta(i)); sin(theta(i))]';
%     s = spSinglePassXYNEW(gt,gammat(i),nbar(i),dnbar(i),rho,2,1);
%     F = s.getAllFish();
%     F1(i) = F(1);
%     F2(i) = F(2);
%     s.getRhoA;
%     s.getRhoAA;
%     MI(i) = s.getMI();
% end