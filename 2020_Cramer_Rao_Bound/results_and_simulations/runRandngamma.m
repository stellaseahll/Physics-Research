clear;
gt = pi/2;
Np = 100;
gammat = logspace(-3,1,50);
nbar = 1;
dnbar = 1e-5*nbar;
Fth = 1./(1+2*nbar).^2./(1+nbar)./nbar;
r = 0.5;
t1 = rand(100,Np)*2*pi;
t2 = rand(100,Np)*2*pi;
p1 = rand(100,Np)*2*pi;
p2 = rand(100,Np)*2*pi;
phi = rand(100,Np)*2*pi;
for j = 1:length(gammat)
    j
    for i = 1:Np
        rho = prepareUncorrelatedState(r,phi(j,i),p1(j,i),t1(j,i),p2(j,i),t2(j,i));
        s = spSinglePassXYNEW(gt,gammat(j),nbar,dnbar,rho,2,1);
        F = real(s.getAllFish());
        F1(j,i) = F(1);
        F2(j,i) = F(2);
        rho = prepareUncorrelatedState(0,phi(j,i),p1(j,i),t1(j,i),p2(j,i),t2(j,i));
        s = spSinglePassXYNEW(gt,gammat(j),nbar,dnbar,rho,2,1);
        F = real(s.getAllFish());
        G1(j,i) = F(1);
        G2(j,i) = F(2);
    end
end