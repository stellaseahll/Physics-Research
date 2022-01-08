function [H F Cint y] = runParametersOLD(Tmin,Tmax,NT,alpha,d)
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
%% F AND H
Fint = dprior.^2./prior;
Fint([1 end]) = 4*pi^2*alpha/(Tmax-Tmin)^2*normfac;
H = trapz(T,prior./T);
F = trapz(T,Fint);

%% HEAT CAPACITY OPTIMIZATION

xguess = 1:Ngap;
options = optimoptions('fmincon','Display','off');
for i = 1:10
    [y(i) fy(i)]= fmincon(@(x) findOptBoundOLD(x,d,prior,T),xguess,[],[],[],[],zeros(1,Ngap)+20*(i-1),zeros(1,Ngap)+20*i,[],options);
   
end

Cint = -fy(fy==min(fy));
y = y(fy==min(fy));
end