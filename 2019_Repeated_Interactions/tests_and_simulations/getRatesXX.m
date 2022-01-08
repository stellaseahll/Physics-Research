% Get effective rates in front of dissipators in analytical ME for pm+mp
% interaction

dimS = 2;
omegaS = 10;
dimB = 2;
omegaB = 10;
delta = omegaS - omegaB;
Omega = omegaS + omegaB;
TB = 10;
gamma = 0.01;
gxx = 1;
gyy = 0;
gzz = 0; 
gpm = 0;
gsb = [gxx gyy gzz gpm];
Ntrials = 100;
rho0 = [1 0; 0 0];

% rhoX = eye(dimS)/dimS;
% rhoTh = exp(-(1:dimS)*omegaB/TB);
% rhoTh = diag(rhoTh)/sum(rhoTh);
% rho2 = exp(-(1:dimS)*omegaS/TB);
% rho2 = diag(rho2)/sum(rho2);
% c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,0.1,gsb,rho0);
% c.findSSwithoutSim();
% 
% H = c.H0{1} + c.Hint{1};
% E = unique(abs(eig(H)));
% rat = E(2)/E(1);
% dtRes = 2*pi/E(1);
% step = (0.01:0.01:4);
dtInt = 1;
% kp = []; %excitation rate (exclude gamma)
% km = []; %deexcitation rate (exclude gamma)
% kz = []; %dephasing rate (exclude gamma)
g = gxx*2;
lambda = sqrt(delta.^2 + g^2);
A = delta + lambda;
N = g^2 + A.^2;
mu = sqrt(Omega.^2 + g^2);
B = Omega + mu;
M = g^2 + B.^2;
p0 = ones(size(omegaB));
p1 = exp(-omegaB/TB);
p0 = p0./(p0+p1);
p1 = 1-p0;
p0 = p0'*ones(size(dtInt));
p1 = p1'*ones(size(dtInt));
k = (2*g*A./N)'*ones(size(dtInt)) .* sin(lambda'*dtInt/2);
m = (2*g*B./M)'*ones(size(dtInt)) .* sin(mu'*dtInt/2);
Num = sum(k.^2.*p0 + m.^2.*p1);
Den = sum(k.^2 + m.^2);
x = Num./Den;