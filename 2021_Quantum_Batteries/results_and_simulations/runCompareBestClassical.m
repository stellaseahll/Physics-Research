%optimal classical parameters
optCtheta = pi/2;
optCq = 0; 
optCE = 1; %units of E: v = sin^2(theta)*(1-2q)

%optimal theoretical quantum
n = 1:3;
for i = 1:length(n)
    [optq(i),y(i)] = fminsearch(@(q) findOptq(q,n(i)),0.25);
end
y = -y;
theta = pi/2./n;

%actual simulation
N = 200;
k0 = 50;
for i = 1:length(n)
%     i
    kmax = k0*n(i);

    [cohE{i},incE{i},~,~] = MERun(optq(i),theta(i),N,kmax);
    
end

for i = 1:length(n)
    kmax = k0*n(i);
    plot(linspace(0,theta(i)*n(i)/(pi/2),kmax+1),cohE{i}(1,:)); hold on;
    plot(linspace(0,theta(i)*n(i)/(pi/2),kmax+1),linspace(0,y(i)*k0,kmax+1)); hold on;
end