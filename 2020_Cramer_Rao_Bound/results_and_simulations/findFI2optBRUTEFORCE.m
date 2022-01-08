function F = findFI2optBRUTEFORCE(theta,gammat,n)

rho = [cos(theta/2) ; sin(theta/2)]*[cos(theta/2) ; sin(theta/2)]';
s = spSinglePassXYZ(pi/2,0,gammat,n,1e-7,rho,2,1); s.alg = 1;
f = s.getAllFish();
F = -real(f(2));