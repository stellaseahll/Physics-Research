function [F Cint y] = runParametersNEW(Tmin,Tmax,NT,alpha,d)
% Tmin = min temp
% Tmax = max temp
% NT = number of T intervals
% alpha = prior parameter
% Np = total number of INTERACTING qubits
% Nlvl = number of distinct energy levels after interaction (excluding
% ground)
Ngap = length(d);
%% TEMPERATURE
dT = (Tmax-Tmin)/NT;
T = Tmin:dT:Tmax; 
%% PRIOR
argm = pi*(T-Tmin)/(Tmax-Tmin);
normfac = 1/(Tmax-Tmin)/(exp(alpha/2)*besseli(0,alpha/2)-1);
prior = normfac*(exp(alpha*(sin(argm)).^2)-1);
dprior = exp(alpha*(sin(argm)).^2)*pi*alpha.*sin(2*argm)/(Tmax-Tmin)*normfac;
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
% %% F
% Fint = dprior.^2./prior;
% Fint([1 end]) = 4*pi^2*alpha/(Tmax-Tmin)^2*normfac;
% F = trapz(T,Fint.*T);

%% HEAT CAPACITY OPTIMIZATION

xguess = 1:Ngap;
options = optimoptions('fmincon','Display','off');
for i = 1:10
[y(i) fy(i)]= fmincon(@(x) findOptBoundNEW(x,d,prior,T),xguess,[],[],[],[],zeros(1,Ngap)+20*(i-1),zeros(1,Ngap)+20*i,[],options);
end
Cint = -fy(fy==min(fy));
y = y(fy==min(fy));
end