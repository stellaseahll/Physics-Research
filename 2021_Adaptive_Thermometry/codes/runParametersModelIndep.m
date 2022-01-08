function [F Cint] = runParametersModelIndep(Tmin,Tmax,NT,alpha,d)
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
M = T.^2.*dprior;
R = myDiff(dT,M);
%% INTERP
Tlong = interp(T,100);
Rlong = interp(R,100);
Mlong = interp(M,100);
idx = find(Rlong<0,1,'last');
Ttilde = Tlong(idx+1);
Cint = abs(Mlong(idx+1));
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
