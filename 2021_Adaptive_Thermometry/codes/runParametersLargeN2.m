function [F const] = runParametersLargeN2(Tmin,Tmax,NT,alpha)

%% TEMPERATURE
dT = (Tmax-Tmin)/NT;
T = Tmin:dT:Tmax; 
%% PRIOR
argm = pi*(T-Tmin)/(Tmax-Tmin);
normfac = 1/(Tmax-Tmin)/(exp(alpha/2)*besseli(0,alpha/2)-1);
prior = normfac*(exp(alpha*(sin(argm)).^2)-1);
dprior = exp(alpha*(sin(argm)).^2)*pi*alpha.*sin(2*argm)/(Tmax-Tmin)*normfac;
M = T.^2.*dprior;
R = myDiff(dT,M);
options = optimoptions('fmincon','Display','off');
minE = 0:100:5000;
minfE = [];
for i = 350:400
i
for j = 1:length(minE)-1
[E(j),fE(j)] =  fmincon(@(x) findOptBoundLargeN(x,i,prior,T),(minE(j+1)-minE(j))/2,[],[],[],[],minE(j),minE(j+1),[],options);
end
minfE = [minfE fE(fE==min(fE))];
end
const = -mean(minfE./(350:400));
%% PRIOR BETA
bmin = 0;
bmax = log(Tmax/Tmin);
db = (bmax-bmin)/NT;
b = bmin:db:bmax;
T2 = exp(b)*Tmin;
argm2 =  pi*(T2-Tmin)/(Tmax-Tmin);
priorb = (exp(alpha*(sin(argm2)).^2)-1);
normfac2 = 1/trapz(b,priorb);
priorb = priorb*normfac2;
dpriorb = myDiff(db,priorb);
Fb = (dpriorb).^2./priorb;
Fb(1) = 4*pi^2*alpha*Tmin^2/(Tmax-Tmin)^2*normfac2;
Fb(end) = 4*pi^2*alpha*Tmax^2/(Tmax-Tmin)^2*normfac2;
F = trapz(b,Fb);