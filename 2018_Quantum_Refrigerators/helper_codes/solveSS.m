
%% Dimensions
for d = 2:4
d
dh = d;
dc = d;
dw = d;

%% Operators for H mode
idh = eye(dh);
Jh = (dh-1)/2;
mh = -Jh:Jh;
t = sqrt((Jh-mh).*(Jh+mh+1));
zh = diag(mh);
ph = diag(t(1:end-1),-1);
mh = ph';

%% Operators for C mode
idc = eye(dc);
Jc = (dc-1)/2;
mc = -Jc:Jc;
t = sqrt((Jc-mc).*(Jc+mc+1));
zc = diag(mc);
pc = diag(t(1:end-1),-1);
mc = pc';

%% Operators for W mode
idw = eye(dw);
Jw = (dw-1)/2;
mw = -Jw:Jw;
t = sqrt((Jw-mw).*(Jw+mw+1));
zw = diag(mw);
pw = diag(t(1:end-1),-1);
mw = pw';

%% Parameters
wh = 5;
wc = 1;
ww = 4;
Th = 2;
Tc = 1;
Tw = 5;
nh = 1/(exp(wh/Th)-1);
nc = 1/(exp(wc/Tc)-1);
nw = 1/(exp(ww/Tw)-1);
g = 0.1;
kh = 0.01;
kc = 0.01;
kw = 0.01;

%% Hamiltonians
Hh = wh*kron(zh,eye(dc*dw));
Hc = wc*kron(eye(dh),kron(zc,eye(dw)));
Hw = ww*kron(eye(dh*dc),zw);
Hint = g*(kron(ph,kron(mc,mw)) + kron(mh,kron(pc,pw)));
H = Hh+Hc+Hw + Hint;

%% Dissipators
Ph = kron(ph,eye(dc*dw));
Pc = kron(eye(dh),kron(pc,eye(dw)));
Pw = kron(eye(dh*dc),pw);
Mh = Ph';
Mc = Pc';
Mw = Pw';
PMh = Ph*Mh;
PMc = Pc*Mc;
PMw = Pw*Mw;
MPh = Mh*Ph;
MPc = Mc*Pc;
MPw = Mw*Pw;
Lh = nh*(lrMultiply(Ph) - 0.5*(leftMultiply(MPh) + rightMultiply(MPh))) + ...
           (nh+1)*(lrMultiply(Mh) - 0.5*(leftMultiply(PMh) + rightMultiply(PMh)));
Lc = nc*(lrMultiply(Pc) - 0.5*(leftMultiply(MPc) + rightMultiply(MPc))) + ...
           (nc+1)*(lrMultiply(Mc) - 0.5*(leftMultiply(PMc) + rightMultiply(PMc)));
Lw = nw*(lrMultiply(Pw) - 0.5*(leftMultiply(MPw) + rightMultiply(MPw))) + ...
           (nw+1)*(lrMultiply(Mw) - 0.5*(leftMultiply(PMw) + rightMultiply(PMw)));

%% Solve Steady
L = 1i*( rightMultiply(Hint)-leftMultiply(Hint)) + kh*Lh + kc*Lc + kw*Lw;
rhoSS = null(L);
rhoSS = reshape(rhoSS,dh*dw*dc,dh*dw*dc);
rho{d} = rhoSS/trace(rhoSS);
rhoSSVec = reshape(rhoSS,(dh*dw*dc)^2,1);
Qh(d) = real(sum(sum(reshape(Lh*rhoSSVec,dh*dw*dc,dh*dw*dc).*Hh)));
Qc(d) = real(sum(sum(reshape(Lc*rhoSSVec,dh*dw*dc,dh*dw*dc).*Hc)));
Qw(d) = real(sum(sum(reshape(Lw*rhoSSVec,dh*dw*dc,dh*dw*dc).*Hw)));
save('data.mat');
end