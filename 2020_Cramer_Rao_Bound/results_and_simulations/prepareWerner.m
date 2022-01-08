function rho = prepareWerner(p)
%entangled for p<0.5
P = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];  %defined for ge basis
Pas = 0.5*(eye(4)-P);
Ps = 0.5*(eye(4)+P);
rho = p* Ps/3 + (1-p)*Pas;
end