% Get effective rates in front of dissipators in analytical ME for pm+mp
% interaction

dimS = 2;
omegaS = 10;
dimB = 2;
omegaB = 0.1:0.1:20;
delta = omegaS - omegaB;
Omega = omegaS + omegaB;
TB = 10;
gamma = 0.01;
gxx = 0;
gyy = 0;
gzz = 0; 
gpm = 0.1;
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
dtInt = 0.01:0.01:10;
% kp = []; %excitation rate (exclude gamma)
% km = []; %deexcitation rate (exclude gamma)
% kz = []; %dephasing rate (exclude gamma)
g = gpm*2;
lambda = sqrt(delta.^2 + g^2);
A = delta + lambda;
N = g^2 + A.^2;
p0 = ones(size(omegaB));
p1 = exp(-omegaB/TB);
p0 = p0./(p0+p1);
p1 = 1-p0;
p0 = p0'*ones(size(dtInt));
p1 = p1'*ones(size(dtInt));
k = (2*g*A./N)'*ones(size(dtInt)) .* sin(lambda'*dtInt/2);
% Num = sum(k.^2.*p0);
% Den = sum(k.^2);
% x = Num./Den
kp = k.^2.*p0;%excitation rate (exclude gamma)
km = k.^2.*p1;%deexcitation rate (exclude gamma)
K = k.^2;
C = g*g*exp(-1i*(delta+lambda)'*dtInt/2)./(N'*ones(size(dtInt))) + exp(-1i*(delta-lambda)'*dtInt/2).*((A.^2.*N)'*ones(size(dtInt)));
kz = abs(1-C).^2; %dephasing rate (exclude gamma)
contourf(omegaS-omegaB,dtInt,log10(K'))
% xlabel('Delta = omegaS-omegaB');
% ylabel('tau');
% %Analytical results
% for j = 1:length(dtInt)
%     k = 2*g*A./N.*sin(lambda*dtInt(j)/2);
%     C = exp(-1i*(delta+lambda)*dtInt(j)/2)*g^2./N + exp(-1i*(delta-lambda)*dtInt(j)/2).*A.^2./N;
%     kz = [kz; abs(1-C).^2];
%     kp = [kp; k.^2.*p0];
%     km = [km; k.^2.*p1];
% end
% 
