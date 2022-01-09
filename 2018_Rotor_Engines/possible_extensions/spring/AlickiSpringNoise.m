function [t meanX meanP meanN meanE] = AlickiSpringNoise(numCyc,Ntrials,nH,nC,w,r,gamma,kappa)
% t = [0:dt:tf]
% gamma = damping rate, kappa = thermalisation rate, w = spring natural
% frequency
% r = effective mass tranlation rate,r=sqrt(hbar*g/(m*L^2)), determines spring stiffness VS pressure

%% Time Intervals
tf = 2*pi/w*numCyc;
dt = tf/numCyc/1000;
t = 0:dt:tf;
Nsteps = length(t);

%% Initialisation
currX = zeros(Ntrials,1);
currP = zeros(Ntrials,1);
currN = zeros(Ntrials,1);
meanX=[0 0];
meanP=[0 0];
meanN=[0 0];
meanE=[0 0];

%% Generation of Random Numbers
dW = normrnd(0,sqrt(dt),Ntrials,Nsteps);

%% Coupling Functions
%valve coupling functions (corresponds to f_H^2 and f_C^2 in our engine, btw 0 and 1)
fH = @(x)( 1 + (abs(x-1)-1).*heaviside(x).*heaviside(2-x) );
fC = @(x)( (abs(x)-1).*heaviside(1-abs(x)) + 1 );

rdt = r^2*dt;
wdt = w^2*dt;
kdt = kappa*dt;
for j = 2:length(t)
    if (mod(j,1000)==1)
        j
    end
    fh = fH(currX);
    fc = fC(currX);
    neff = (fh * nH + fc * nC)./(fh + fc);
    keff = kappa*(fh+fc);
    tmpX = currX + currP*dt;
    tmpP = currP + (currN).*rdt - currX.*wdt - gamma.*currP.*dt; %+ backaction noise  ;
    tmpN = currN + ((nH-currN).*fh + (nC-currN).*fc)*kdt + sqrt(2*currN.*keff.*neff).*dW(:,j); %checked, noise terms same as rotor engine
    tmpN(tmpN<0) = 0; %change to N = a + ib would be better
    currX = tmpX;
    currP = tmpP;
    currN = tmpN;
    meanN(j,:) = [mean(currN) std(currN)];
    meanX(j,:) = [mean(currX) std(currX)];
    meanP(j,:) = [mean(currP) std(currP)];
    meanE(j,:) = [mean(currX.^2 + currP.^2) std(currX.^2 + currP.^2)];

end