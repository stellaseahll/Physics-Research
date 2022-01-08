function F = findFI3optBRUTEFORCE(theta,gammat,n)

rho = [cos(theta/2) ; sin(theta/2)]*[cos(theta/2) ; sin(theta/2)]';
s = spSinglePassXYZ(pi/2,0,gammat,n,1e-7,rho,3,1); s.alg = 2;
f = s.getAllFish();
F = -real(f(3));