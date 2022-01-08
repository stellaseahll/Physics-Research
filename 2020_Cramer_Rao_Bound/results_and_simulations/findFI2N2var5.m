function F = findFI2N2var5(x,gammat,n)

theta1 = x(1);
theta2 = x(2);
phi = x(3);
alpha = x(4);
r = abs(cos(x(5)));
ket1p = [cos(theta1/2); sin(theta1/2)];
ket2p = [cos(theta2/2); exp(1i*phi)*sin(theta2/2)];
ket1m = [sin(theta1/2); -cos(theta1/2)];
ket2m = [exp(-1i*phi)*sin(theta2/2); -cos(theta2/2)];
ket = sqrt(r)*kron(ket1p,ket2p)+sqrt(1-r)*kron(ket1m,ket2m)*exp(1i*alpha);
rho = ket*ket';
s = spSinglePassXYZ(pi/2,0,gammat,n,1e-6*n,rho,2,1); s.alg = 1;
f = s.getAllFish();
F = -real(f(2)); %INSERT YOUR FUNCTION TO GET F2 HERE, BUT TAKE NEGATIVE VALUE.