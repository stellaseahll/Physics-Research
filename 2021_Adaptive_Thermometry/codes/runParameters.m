function runParameters(Tmin,Tmax,NT,alpha,N,Np,Nlvl)
% Tmin = min temp
% Tmax = max temp
% NT = number of T intervals
% alpha = prior parameter
% N = total number of qubits
% Np = total number of INTERACTING qubits
% Nlvl = number of distinct energy levels after interaction 

dT = (Tmax-Tmin)/NT;



Nset = ceil(N/Np);
Eguess = [2:Nlvl ones(1,Nlvl-2) 2^Np-Nlvl+1]; %non-degenerate spectrum
options = optimoptions('fmincon','Display','off');

for i = 1:100
    if (mod(i,100) ==1)
        i
    end
%     nonlcon = @constraints;
    [y(i,:) fy(i)]= fmincon(@(x) findOptBoundwDegeneracy(x,Tmin,Tmax,NT,i,Nlvl,Np,alpha),Eguess,[],[],[zeros(1,Nlvl-1) ones(1,Nlvl-1)],2^Np-1,[zeros(1,Nlvl-1) ones(1,Nlvl-1)],[Inf*ones(1,Nlvl-1) 2^Np-1*ones(1,Nlvl-1)],[],options);
    E(i,:) = [0 y(i,1:Nlvl-1)];
    d(i,:) = round([1 y(i,Nlvl:end)]);
    d(i,end) = 2^Np - sum(d(i,1:end-1));
end