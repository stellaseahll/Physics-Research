function [t W Qh Qc] = getTimeEvolve(wp)

%Analytical solution for:
%fh = [0 0.75pi], fc = [pi 1.75pi], Hon = [0 pi]

numCyc = 200;
sz = -0.5;
Qh = 0; 
Qc = 0; 
W = 0;
% wp = 1;
nh = 1; 
nc = 0.01;
kh = 0.01; 
kc = 0.01;
g = 0.01;
ws = 1;

t = 0;
thetaH0 = 0;
thetaH1 = 0.75*pi;
thetaH = linspace(thetaH0,thetaH1,51);
thetaH(1) = [];
thetaC0 = pi;
thetaC1 = 1.75*pi;
thetaC = linspace(thetaC0,thetaC1,51);
thetaC(1) = [];

for i = 1:numCyc
    p0 = sz(end);
    E0 = (ws+g)*sz(end);
    %from 0 to thetaH0
    t(end+1) = (i-1)*2*pi + thetaH0;
    sz(end+1) = p0;
    Qh(end+1) = Qh(end);
    Qc(end+1) = Qc(end);
    W(end+1) = W(end) + p0*(-g);
    %from thetaH0 to thetaH1
    t = [t (i-1)*2*pi+thetaH];
    ph = ((nh+1)-(2*nh+1)*(0.5-p0))*exp(-kh*(2*nh+1)*(thetaH-thetaH0)/wp)/(2*nh+1)-0.5/(2*nh+1);
    qh = (ws+g)*(((nh+1)-(2*nh+1)*(0.5-p0))*exp(-kh*(2*nh+1)*(thetaH-thetaH0)/wp)/(2*nh+1)) + Qh(end) -...
        (ws+g)*(((nh+1)-(2*nh+1)*(0.5-p0))/(2*nh+1));
    sz = [sz ph];
    Qh = [Qh qh];
    Qc = [Qc Qc(end)+zeros(size(ph))];
    W = [W W(end)+zeros(size(ph))];
    %from thetaH1 to thetaC0
    t(end+1) = (i-1)*2*pi + thetaC0;
    p0 = sz(end);
    sz(end+1) = p0;
    Qh(end+1) = Qh(end);
    Qc(end+1) = Qc(end);
    W(end+1) = W(end) + p0*(g);
    %from thetaC0 to thetaC1
    t = [t (i-1)*2*pi+thetaC];
    pc = ((nc+1)-(2*nc+1)*(0.5-p0))*exp(-kc*(2*nc+1)*(thetaC-thetaC0)/wp)/(2*nc+1)-0.5/(2*nc+1);
    qc = (ws)*(((nc+1)-(2*nc+1)*(0.5-p0))*exp(-kh*(2*nc+1)*(thetaC-thetaC0)/wp)/(2*nc+1)) + Qc(end) -...
        (ws)*(((nc+1)-(2*nc+1)*(0.5-p0))/(2*nc+1));
    sz = [sz pc];
    Qc = [Qc qc];
    Qh = [Qh Qh(end)+zeros(size(pc))];
    W = [W W(end)+zeros(size(pc))];
end    
    