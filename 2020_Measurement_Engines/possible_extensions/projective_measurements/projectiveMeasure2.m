% 26/09/2018 Projective Measurement

%% Define Initial state and basis to be projected
clear;clc;clf;
% clf;plot(cos(0:0.01:2*pi),sin(0:0.01:2*pi),'k--'); hold on; axis equal;
% Pauli matrices
sz = [1 0; 0 -1];
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
seye = [1 0; 0 1];

Theta = (-0.5:0.01:0.5)*pi;
c = colormap(lines); R = 1;
Alpha = 2;%logspace(-2,0,100);
for alpha = 1:length(Alpha)
for i = 1:length(Theta)
% Initial state vector
r = 1;
theta = Theta(i);
phi = 0;
rho1 = createState([r,theta,phi]);
[bvec1, p] = eig(rho1);  %bvec = basis vectors of rho, eig = corresponding eigenvalues
p = diag(p);
bvec1 = bvec1(:,end:-1:1); %arranged in descending order
p = p(end:-1:1); %arranged in descending order
U1 = bvec1'; % U = to rotate from energy basis to bvec

% Basis of projection
Ptheta = pi/2;
Pphi = 0;
[bvec2, eig2] = eig(createState([1,Ptheta,Pphi]));  %bvec = basis vectors of rho, eig = corresponding eigenvalues
U2 = bvec2'; % U = to rotate from energy basis to bvec
rho4 = zeros(2);
P = zeros(2,1);
for k = 1:2
    P(k) = 0;
    for j = 1:2
        P(k) = P(k) + (p(j) * abs(bvec2(:,k)' * bvec1(:,j))^2); 
    end
    rho4 = rho4 +  P(k)*(bvec2(:,k)*bvec2(:,k)');
end
[P,idx] = sort(P,'descend');
bvec2 = bvec2(:,idx);
% 
% plotState(rho1);
% plotState(createState([1,Ptheta,Pphi]))
% plotState(rho4);
% Hamiltonian 
H = -sz;

%% Compute work and heat
E = -log(P);
idx = find(P>eps*1e5);
idx2 = find(p>eps*1e5);
Q(i) = sum((P(idx)-p(idx)).*E(idx));
U(i) = trace(H*(rho4-rho1))*Alpha(alpha);
W(i) = U(i)-Q(i);
Q2(i) = -(sum(P(idx).*log(P(idx))) - sum(p(idx2).*log(p(idx2))));
W2(i) = U(i)-Q2(i);
end
Rat(alpha) = trapz(U)/trapz(Q);
cumW(alpha) = trapz(Theta,W);
plot(Theta/pi,U,'color',c(R,:)); hold on;
plot(Theta/pi,Q2,'color',c(R,:),'linestyle',':'); hold on;
plot(Theta/pi,W2,'color',c(R,:),'linestyle','--'); hold on;
% % plot(Theta/pi,real(W)./real(Q),'color',c(R,:)); hold on;
% R = R + 1;
end

% plot(Theta/pi,real(Q),'k--');
