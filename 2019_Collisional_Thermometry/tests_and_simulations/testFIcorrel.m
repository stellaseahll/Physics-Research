%random parameter draws, see correlation plots

N = 2000;
nbar = 1+ rand(N,1)*0;
%gt = pi/2*10.^(-2+2*rand(N,1));
gt = 0*pi/100 + pi/2*rand(N,1);
gammat = 0.5 + 0*rand(N,1) + 0*10.^(-2+2*rand(N,1));
theta = linspace(0,2,N).'+0*rand(N,1); %state Bloch angle in units pi.
p0 = 1 + 0*linspace(0,1,N).' + 0.0*rand(N,1);

dnbar=0.0001;
np=3;

F = zeros(N,np); %Q FI
f = zeros(N,np); %class FI
dpA = zeros(N,1); %change in probe pop
C = zeros(N,1); %correl between subsequent probes
I = zeros(N,1); %mutual info
Dg = zeros(N,1); %mutual info
D = zeros(N,1); %mutual info
dC = zeros(N,1); %diff by nbar
ddpA = zeros(N,1); 
tic
for n=1:N
    %rhop = [1+cos(theta(n)*pi) sin(theta(n)*pi); sin(theta(n)*pi) 1-cos(theta(n)*pi)]/2;
    rhop = [p0(n), 0; 0, 1-p0(n)];
    s = spSinglePassXY(gt(n),gammat(n),nbar(n),dnbar,rhop,np,1);
    s.alg=3;
    [FF,ff] = s.getAllFish();
    M = s.getAllSignals(2);
    M1 = s.getAllSignals(1);
    M3 = s.getAllSignals(3);
    s.getRhoAA();
    s.getRhoA();
    F(n,:)=FF;
    f(n,:)=ff;
    dpA(n) = mean(diag(M));
    C(n) = mean(diag(M,-1));
    dC(n) = mean(diag(M3,-1)-diag(M1,-1))/2/dnbar;
    ddpA(n) = mean(diag(M3)-diag(M1))/2/dnbar;
    I(n) = 2*s.getSvn(s.rhoA{2})-s.getSvn(s.rhoAA{2});
    Dg(n) = s.getGeometricDiscord();
    D(n) = s.getDiscordEnergyBasis();
end
toc

%plot some data points to see link between correlations and FI
%enhancement...
%figure(); loglog(abs(C./dpA),real(F(:,end)-f(:,end))./real(F(:,end)+f(:,end)),'.')

%plot(abs(I),real(F(:,2)-f(:,2))./real(F(:,2)+f(:,2)),'.',...
%    abs(I),real(D)*100,'.',abs(I),abs(C)/max(C),'.')

%plot(D,real(F(:,3)-f(:,3))./real(F(:,3)+f(:,3)),'.',...
%    D,real(I),'.',D,abs(C)/max(C),'.',D,F(:,3)./F(:,2)*2/3,'.')

%plot(Dg,real(F(:,3)-f(:,3))./real(F(:,3)+f(:,3)),'.',...
%    Dg,real(I),'.',Dg,abs(C)/max(C),'.',Dg,F(:,3)./F(:,2)*2/3,'.')
