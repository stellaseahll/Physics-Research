
%% Qubit Operators
sz = 0.5*[-1 0 ; 0 1];
sp = [0 0; 1 0];
sm = [0 1; 0 0];
sx = sp+sm;
id = eye(2);
sz1 = kron(sz,id);
sz2 = kron(id,sz);
%% Rates/Timescales
w1 = 10; % spin 1 freq
w2 = 10; %spin 2 freq
g = 1; % Coupling between system spins
k = 1; % Coupling between bath+system spin
gamma = 1e-3; % Frequency of interaction
tInt = 0.00001/gamma; %interaction time tInt << gamma^-1
tStep = 0.5/gamma; % tInt << tStep << gamma^-1
nEnd = 50000; %relative end time
tEnd = nEnd/gamma;
Nsteps = nEnd*2; %number of steps to get to tEnd

%% Other parameters
T1 = 10; %bath 1 temp
T2 = 10; %bath 2 temp
wb1 = w1; %bath 1 freq
wb2 = w2; %bath 2 freq
Ntrials = 1; 
%% Draw wait time for jump
% Distribution: gamma*exp(-gamma*tWait), tWait: waiting time
% Cumulative distribution: 1- exp(-gamma*tWait)
p = rand(Ntrials,Nsteps); 
tWait = log(1-p)/(-gamma);
tWait = cumsum(tWait,2);
tWait = tWait +  (ones(Ntrials,1)*(0:Nsteps-1))*tInt;

%% Hamiltonians
Hfree = w1*kron(sz,id) + w2*kron(id,sz) + g*kron(sx,sx);%
[v1 e1] = eig(Hfree);
rhoG = zeros(4);
for n = 1:4
        rhoG = rhoG + v1(:,n)*v1(:,n)'*exp(-e1(n,n)/T1);
end
rhoG = rhoG/trace(rhoG);
Ufree = zeros(4);
for n = 1:4
    V{n} =  v1(:,n)*v1(:,n)';
    Ufree = Ufree + V{n}*exp(-1i*e1(n,n)*tStep); %Ufree = free evolution for tStep
end
Htotal = wb1*kron(id,kron(id,sz))+ kron(Hfree,id)+k*kron(sx,kron(id,sx));
[v2 e2] = eig(Htotal);
Uint = zeros(8);
for j = 1:8
    W{j} = v2(:,j)*v2(:,j)';
    Uint = Uint + W{j}*exp(-1i*e2(j,j)*tInt);
end

%% States
rhoS1 = [1 0; 0 0];
rhoS2 = [1 0; 0 0];
rhoS = kron(rhoS1,rhoS2);
rhoB1 = [1 0; 0 exp(-wb1/T1)];
rhoB1 = rhoB1/trace(rhoB1);
rhoB2 = [1 0; 0 exp(-wb2/T2)];
rhoB2 = rhoB2/trace(rhoB2);
rhoP = kron(rhoB1,rhoB2);
z1 = zeros(Ntrials,Nsteps);
z2 = z1;
t = 0:tStep:tEnd;
for j = 1:Ntrials
    idx = 1;
    z1(j,idx) = trace(rhoS*sz1);
    z2(j,idx) = trace(rhoS*sz2);
    for m = 1:Nsteps
        if (idx==length(t))
            break;
        end
        for n = 1:floor(tWait(j,m)/tStep)
            idx = idx + 1;
            rhoS = Ufree*rhoS*Ufree';
            z1(j,idx) = trace(rhoS*sz1);
            z2(j,idx) = trace(rhoS*sz2);
            if (idx==length(t))
                break;
            end
        end
        if (idx==length(t))
           break;
        end
        idx = idx + 1;

        tbefore = mod(tWait(j,m),tStep);
        tafter = ceil(tWait(j,m)/tStep)*tStep - tWait(j,m) - tInt;
        if (tafter>=0)
            Ubefore = zeros(4);
            Uafter = zeros(4);
            for n = 1:4
                Ubefore = Ubefore + V{n}*exp(-1i*e1(n,n)*tbefore);
                Uafter = Uafter + V{n}*exp(-1i*e1(n,n)*tafter);
            end
            rhoS = Uafter*ptrace(Uint*kron(Ubefore*rhoS*Ubefore',rhoB1)*Uint',4,2)*Uafter';
            z1(j,idx) = trace(rhoS*sz1);
            z2(j,idx) = trace(rhoS*sz2);
            tWait(j,m+1:end) = tWait(j,m+1:end)-ceil(tWait(j,m)/tStep)*tStep;
        else
            Ubefore = zeros(4);
            for n = 1:4
                Ubefore = Ubefore + V{n}*exp(-1i*e1(n,n)*tbefore);
            end
            tInt1 = (tStep-tbefore);
            tInt2 = tInt - tInt1;
            Uint1 = zeros(8);
            for n = 1:8
                Uint1 = Uint1 + W{n}*exp(-1i*e2(j,j)*tInt1);
            end
            rhotmp = Ubefore*rhoS*Ubefore';
            rhoS = ptrace(Uint1*kron(rhotmp,rhoB1)*Uint1',4,2);
            z1(j,idx) = trace(rhoS*sz1);
            z2(j,idx) = trace(rhoS*sz2);
            rhoS = ptrace(Uint*kron(rhotmp,rhoB1)*Uint',4,2);
            tWait(j,m+1:end) = tWait(j,m+1:end)-ceil(tWait(j,m)/tStep)*tStep-tInt2;
        end
    end
end
clf;
plot(t*gamma,real(z1),'b-',t*gamma,real(z2),'r-'); hold on;
plot([t(1) t(end)]*gamma,[trace(rhoB1*sz) trace(rhoB1*sz)],'k--');
hold off;
rhoP
rhoG
rhoS
