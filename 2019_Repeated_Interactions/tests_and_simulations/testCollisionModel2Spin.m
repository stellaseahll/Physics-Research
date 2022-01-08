%test script single spin coupled to bath
%Load the object with inputs
%Arguments:
% (dimS,omegaS,dimB,omegaB,TB,gamma,dtint,gxx,gyy,gzz,gpm)
%  -> gamma = overall jump rate
%  -> dtint = length of one collision
%  -> gxx,gyy,gzz,gpm coupling rates in interaction Hamiltonian:
%       Hint = gxx*JxS*JxB + gyy*JyS*JyB + gzz*JzS*JzB + gpm*(J-S*J+B + h.c.)
clear;clc;
% Parameters
dimS =3;
omegaS = 10;
omegaB = 8;
dimB = 2;
J = (dimS-1)/2;
TB = 10;
gxx1 = 1;
gyy1 = 1;
gzz1 = 0; 
gpm1 = 0;
gxx2 = 0;
gyy2 = 0;
gzz2 = 0; 
gpm2 = 0;
gxx = 1;
gyy = 1;
gzz = 0; 
gpm = 0;
gsb1 = [gxx1 gyy1 gzz1 gpm1];
gsb2 = [gxx2 gyy2 gzz2 gpm2];
gss = [gxx gyy gzz gpm];
dtInt = logspace(0,2,200);
rho0 = zeros(dimS*dimS);
rho0(1) = 1;
sz = diag(-J:J);
% sp = kron([0 0; 1 0],eye(dimS));
% sm = kron([0 1; 0 0],eye(dimS));
obs{1} = kron(sz,eye(dimS));
obs{2} = kron(eye(dimS),sz);
obs2{1} = sz;
gamma = 1e-3;
dtStep = 10/gamma;
tEnd = 1000*dtStep;
t = 0:dtStep:tEnd;
% c1 = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt,gsb1);
% [rhoSt1, Jzt1, NjumpsAv1, t1] = c1.runTheBitch(2000,50,200);
% c = collisionModel2SpinBath(dimS,omegaS,dimS,omegaS,dimB*[1 1],omegaB*[1 1],TB*[1 1],gamma*[1 0],dtInt(1)*[1 1],gsb1,gsb2,gss,rho0);
% for i = 1:length(obs)
%     obs{i} = c.basisTransform(obs{i});
% end
% HS = c.HS;
% eigHS = eig(HS);
% emax = max(eigHS);
% emin = abs(min(eigHS));
% omegaB = 10;%[emax-emin, emax+emin, 2*emax, 2*emin];
% rhoth = expm(-HS/TB);
% rhoth = expm(-HS/TB)/trace(rhoth);
rhoth1 = exp(-omegaB*(-J:J)/TB);
rhoth1 = diag(rhoth1/sum(rhoth1));
rhoth1 = kron(rhoth1,rhoth1);
% rhomix = kron(eye(dimS)/dimS,eye(dimS)/dimS);

for j = 1:length(dtInt)
%         L = c.getHdiss();
%         L2 = zeros(dimS^2);
%         for k = 1:length(omegaB)
            c = collisionModel2SpinBath(dimS,omegaS,dimS,omegaS,dimB*[1 1],omegaB*[1 1],TB*[1 1],gamma*[1 0],dtInt(j)*[1 1],gsb1,gsb2,gss,rho0);
            c.findSSwithoutSim();
           
            rho{j} = c.Ucomp2E'*c.rhoSSscatter{1}*c.Ucomp2E;
            z1(j) = sum(sum(rho{j}.*obs{1}));
            z2(j) = sum(sum(rho{j}.*obs{2}));
%         distTh(j) = norm(rho{j}{1}-rhoth);
%         [a b] = ptrace(rho{j}{1},dimS,dimS);
%         dist2Spin(j) = norm(a-b);
%         distProd(j) = norm(rho{j}{1}-kron(a,b));
            A = (rho{j}-rhoth1);
            distProdTh(j) = sqrt(trace(A'*A));
%         distMix(j) = norm(rho{j}{1}-rhomix);
% 	    [rhoStME{j}, obstME{j}] = c.runMEScatter(t,obs);
%         c2 = collisionModelSpinBath(dimS,omegaS,dimB,omegaB,TB,gamma,dtInt(i),diag([gsb1 0]),[1 0; 0 0]);
%         c2.prepareSim;
%         [rhoSt2ME{j}, obstME2{j}] = c2.runMEScatter(t,obs2);
%         L = -1i*(c.leftMultiply(c.HS)-c.rightMultiply(c.HS)) + gamma*c.pthB{1}(1)*(c.rightMultiply(sp)*c.leftMultiply(sm) - 0.5*c.leftMultiply(sp*sm) - 0.5*c.rightMultiply(sp*sm)) +...
% gamma*c.pthB{1}(2)*(c.rightMultiply(sm)*c.leftMultiply(sp) - 0.5*c.leftMultiply(sm*sp) - 0.5*c.rightMultiply(sm*sp));
%         localRho{j}(:,:,1) = rho0;
%         Zt1(1) = sum(sum(localRho{j}(:,:,1).*obs{1}));
%         Zt2(1) = sum(sum(localRho{j}(:,:,1).*obs{2}));
%         for k = 2:length(t)
%             localRho{j}(:,:,k) = reshape(expm(L*t(k))* reshape(rho0,dimS^4,1),dimS^2,dimS^2);
%             Zt1(k) = trace(localRho{j}(:,:,k)*kron(sz,eye(2)));
%             Zt2(k) = trace(localRho{j}(:,:,k)*kron(eye(2),sz));
%         end
end
% semilogx(dtInt,z1,'k',dtInt,z2,'ko');
clf;
semilogx(dtInt,distProdTh);
% figure;
% semilogx(dtInt,z1,dtInt,z2);
% plot(t(1:20:end),obstME2{1}(1,1:20:end),'r^-')
% hold on;
% p0 = c.pthB{1}(1);
% Jzth = p0*-0.5 + (1-p0)*0.5;
% plot([0 t(end)], [Jzth Jzth],'--');