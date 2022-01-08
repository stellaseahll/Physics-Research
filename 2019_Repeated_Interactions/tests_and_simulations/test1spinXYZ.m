%test script single spin coupled to bath with analytical expression for
%general linear coupling
%Load the object with inputs
% Parameters
dimS = 2;
omegaS = 10;
dimB = 2;
delta = 0;
omegaB = omegaS-delta;
Omega = omegaS + omegaB;
TB = 10;
gamma = 1e-3;
gxx = 1;
gyy = 0.4;
gzz = 0;
g(1,1) = gxx/2;
g(2,2) = gyy/2;
g(3,3) = gzz/2;
g(5,5) = 0;
Ntrials = 100;
rho0 = [1 0; 0 0];
dtInt = 0.01;
% rhoSS = zeros(dimS,dimS,length(dtInt));
% t = 0:100000:1000000;
% obs{1} = [-1 0; 0 1];
% obs{2} = [0 1; 1 0];
% for j = 1:length(dtInt)
%     c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(j),g,rho0);
%     c.findSSwithoutSim;
%     [rhoSt{j}, obst{j}] = c.runMEScatter(t,obs);
% end
% Analytical Results
Gm = (gxx - gyy)/4;
Gl = (gxx + gyy)/4;
mSq = sqrt(Omega^2+Gm^2);
mP = gzz/4+mSq;
mM = gzz/4-mSq
mTmp = Omega + mSq;
mNorm = (Gm)^2 + mTmp^2;
lSq = sqrt(delta^2+Gl^2);
lP = -gzz/4+lSq;
lM = -gzz/4-lSq;
lTmp = delta + lSq;
lNorm = (Gl)^2 + lTmp^2;
C = exp(-1i*delta*dtInt/2)*(exp(-1i*lSq*dtInt/2)*Gl^2 + exp(1i*lSq*dtInt/2)*lTmp^2)/lNorm;
D = exp(-1i*Omega*dtInt/2)*(exp(-1i*mSq*dtInt/2)*Gm^2 + exp(1i*mSq*dtInt/2)*mTmp^2)/mNorm;
m = 2*Gm*mTmp*sin(mSq*dtInt/2)*exp(-1i*Omega*dtInt/2)/mNorm;
k = 2*Gl*lTmp*sin(lSq*dtInt/2)*exp(-1i*delta*dtInt/2)/lNorm;
theta = angle(k*m);
R = abs(k*m);
sx = [0 0.5; 0.5 0];
sy = [0 1i/2; -1i/2 0];
sz = [-0.5 0; 0 0.5];
sp = [0 0; 1 0];
sm = [0 1; 0 0];
s0 = [1 0; 0 0];
s1 = [0 0; 0 1]; 
S = (sx*cos(theta/2)+sy*sin(theta/2));
rhoTh = exp(-(1:dimS)*omegaB/TB);
rhoTh = diag(rhoTh)/sum(rhoTh);
p0 = rhoTh(1);
p1 = rhoTh(4);
for j = 1:length(dtInt)
    c = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(j),g,rho0);
    c.findSSwithoutSim();
    c.getLpm;    
%     rhoS(:,:,j) = c.rhoSS{1}; 
%     z1(j) = rhoS(1,1,j);
%     rhoSscatter(:,:,j) = c.rhoSSscatter{1};
%     z2(j) = rhoSscatter(2,2,j);
%     c.getSteadyVal();
%     trDistX(j)  = norm(rhoS(:,:,j)-rhoX);
%     trDistTh(j) = norm(rhoS(:,:,j)-rhoTh);
%     trDist(j) = norm(rhoS(:,:,j)-rhoSscatter(:,:,j));
    K0 = -1i*omegaS*(c.leftMultiply(sz)-c.rightMultiply(sz)) ;
    K1 = (k^2*p0+m^2*p1-R)*(c.rightMultiply(sp)*c.leftMultiply(sm) - 0.5*c.leftMultiply(sp*sm) - 0.5*c.rightMultiply(sp*sm)) +...
        (k^2*p1+m^2*p0-R)*(c.rightMultiply(sm)*c.leftMultiply(sp) - 0.5*c.leftMultiply(sm*sp) - 0.5*c.rightMultiply(sm*sp));
    K2 = abs(D-C)^2*(c.rightMultiply(sz)*c.leftMultiply(sz) - 0.25*eye(4)) + 4*R*(c.rightMultiply(S)*c.leftMultiply(S) - 0.5*c.rightMultiply(S'*S) -0.5*c.leftMultiply(S'*S));
    K3 = 1i*imag(C*D)*(c.leftMultiply(sz)-c.rightMultiply(sz));
    K = K0 + gamma*(K1+K2-K3);
    rhoCal = null(K);
    rhoCal = reshape(rhoCal,dimS,dimS)/sum(rhoCal)
    c.rhoSSscatter{1}
end