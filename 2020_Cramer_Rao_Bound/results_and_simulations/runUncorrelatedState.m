clear;
gt = pi/2;
for nbar = [1 5 10 15 20]
% Gamma = 1;
gammat = 0.1;

dnbar = 1e-5;
Fth = 1/(1+2*nbar)^2/(1+nbar)/nbar;
N  = 0:0.05:0.5;
r = (1-sqrt(1-2*N))/2;
Np = 10000;
t1 = rand(length(r),Np)*2*pi;
t2 = rand(length(r),Np)*2*pi;
p1 = rand(length(r),Np)*2*pi;
p2 = rand(length(r),Np)*2*pi;
phi = rand(length(r),Np)*2*pi;
    for k = 1:length(r)
        for i = 1:Np
            if (mod(i,1000) == 0)
                fprintf('(%d,%d)\n',k,i);
            end
            rho = prepareUncorrelatedState(r(k),phi(k,i),p1(k,i),t1(k,i),p2(k,i),t2(k,i));
            s = spSinglePassXYNEW(gt,gammat,nbar,dnbar,rho,2,1);
            F = real(s.getAllFish());
            F1(k,i) = F(1);
            F2(k,i) = F(2);
        end
    end
filename = sprintf('data_nbar%d_gammat%.2f.mat',nbar,gammat);
save(filename);
end