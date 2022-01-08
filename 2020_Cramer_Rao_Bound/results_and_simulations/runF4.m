clear;clc;
load('F2N2.mat');
[maxF2 gammaidx] = max(b2_F2);
Nt = 5;
% N4_F4 = 2*maxF2;

for i = 1:length(gammat)
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        x0 = [b2_theta1(gammaidx(i),i); b2_theta2(gammaidx(i),i);b2_phi(gammaidx(i),i);b2_alpha(gammaidx(i),i);acos(b2_r(gammaidx(i),i))];
        N4_F4(i,j) = -findFI2N4var5(x0,gammat(i),nbar(j));
    end
end
% for i = 1:length(nbar)
%     
%     options = optimset('MaxFunEvals',1e5,'TolFun',10^min([floor(log10(2*maxF2(i)))-4,-6]));
%     phi = b2_phi(gammaidx(i),i)*(1+randn(1,Nt)*0.05);
%     alpha = b2_alpha(gammaidx(i),i)*(1+randn(1,Nt)*0.05);
%     theta1 = b2_theta1(gammaidx(i),i)*(1+randn(1,Nt)*0.05);
%     theta2 = b2_theta2(gammaidx(i),i)*(1+randn(1,Nt)*0.05);
%     coeff = acos(b2_r(gammaidx(i),i))*(1+randn(1,Nt)*0.05);
%     N4_phi(i) = b2_phi(gammaidx(i),i);
%     N4_alpha(i) = b2_alpha(gammaidx(i),i);
%     N4_theta1(i) = b2_theta1(gammaidx(i),i);
%     N4_theta2(i) = b2_theta2(gammaidx(i),i);
%     N4_coeff(i) = acos(b2_r(gammaidx(i),i));
%     x0 = [b2_theta1(gammaidx(i),i); b2_theta2(gammaidx(i),i);b2_phi(gammaidx(i),i);b2_alpha(gammaidx(i),i);acos(b2_r(gammaidx(i),i))];
%     N4_F4(i) = -findFI2N4var5(x0,gammat(gammaidx(i)),nbar(i));
% %     for j = 1:Nt
% %         fprintf('(%d,%d)\n',i,j);
% %         xguess = [theta1;theta2;phi;alpha;coeff];
% %         [x1, f1] = fminsearch(@(x) findFI2N4var5(x,gammat(gammaidx(i)),nbar(i)),xguess,options);
% %         Ftmp = -f1;
% %         if (Ftmp>N4_F4(i))
% %             N4_F4(i) = Ftmp;
% %             N4_phi(i) = x1(1);
% %             N4_alpha(i) = x1(2);
% %             N4_theta1(i) = x1(3);
% %             N4_theta2(i) = x1(4);
% %             N4_coeff(i) = x1(5);
% %         end
% %     end
% end

