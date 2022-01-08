sz = [1 0; 0 -1];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
seye = [1 0; 0 1];
Theta = pi/4;
for k = 1:length(Theta)
% Initial state vector
r = 1;
theta = Theta(k);
phi = pi/2;
ny = r*sin(theta)*sin(phi);
nx = r*sin(theta)*cos(phi);
nz = r*cos(theta);

rhoIn = (seye + nx*sx + ny*sy + nz*sz)/2;

% Axis to be projected on
Pr = 1;
Ptheta = pi/2;
Pphi = 0;
Py = Pr*sin(Ptheta)*sin(Pphi);
Px = Pr*sin(Ptheta)*cos(Pphi);
Pz = Pr*cos(Ptheta);


rhoP =  (seye + Px*sx + Py*sy + Pz*sz)/2;
[Pv Peig] = eig(rhoP);
U = zeros(2,2);

for j = 1:2
    U = U + seye(:,j)*Pv(:,j)';
end

% Rotate state
rhoRot = U*rhoIn*U';
mz = (rhoRot(4) - rhoRot(1));
mx = real(rhoRot(1,2))*2;
my = imag(rhoRot(1,2))*2;

% Dephase
rhoDep = diag(diag(rhoRot));
dz = (rhoDep(4) - rhoDep(1));
dx = real(rhoDep(1,2))*2;
dy = imag(rhoDep(1,2))*2;

% Rotate back
rhoF = U'*rhoDep*U;
fz = (rhoF(4) - rhoF(1));
fx = real(rhoF(1,2))*2;
fy = imag(rhoF(1,2))*2;

W(k) = trace(rhoF*sz) - trace(rhoIn*sz);
pF = abs(eig(rhoF));
pI = abs(eig(rhoIn));
pI(pI<=10*eps) = [];
pF(pF<=10*eps) = [];
S(k) = sum(-pF.*log(pF)) + sum(pI.*log(pI));
Q(k) = 1-trace(rhoIn*sz).^2;
end
plot(W); hold on;
plot(S);
plot(Q); hold off;
clf;
hold on; axis equal;
% % [X, Y, Z] = sphere;
% % plot3(X,Y,Z);
plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k--');
plot([0 nx],[0,nz],'b','linewidth',2);
plot([0 Px],[0,Pz],'r','linewidth',2);
plot([0 fx],[0,fz],'g--','linewidth',2);
