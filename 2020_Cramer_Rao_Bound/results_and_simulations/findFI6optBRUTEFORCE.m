function F = findFI6optBRUTEFORCE(theta,gammat,n)

rho = [cos(theta/2) ; sin(theta/2)]*[cos(theta/2) ; sin(theta/2)]';
s = spSinglePassXYZ(pi/2,0,gammat,n,1e-5,rho,6,1); s.alg =2r;
f = s.getAllFish();
F = -real(f(6));