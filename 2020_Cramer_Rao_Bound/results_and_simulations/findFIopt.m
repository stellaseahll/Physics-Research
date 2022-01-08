function F = findFIopt(theta,gammat,n)
theta = pi-theta;
L = 1-exp(-gammat*(2*n+1));
rhoA(1,1) = -((-1 - 2*n + L + (1 + 2*n)*(-1 + L)*cos(theta))/(2 + 4*n));
rhoA(2,2) = 1-rhoA(1,1);
rhoA(1,2) = (sqrt(1 - L)*cos(theta)*(L + (1 + 2*n)*(-1 + L)*cos(theta))*sin(theta))/(2 + 4*n);
rhoA(2,1) = rhoA(1,2);
drhoA(1,1) = (exp(-(1 + 2*n)*gammat)*(-1 + exp(gammat + 2 *n *gammat) ...
    - gammat -  2* n* gammat - (1 + 2* n)^2* gammat*cos(theta)))/(1 + 2*n)^2;
drhoA(2,2) = -drhoA(1,1);
drhoA(1,2) = ((exp(-(1 + 2*n)*gammat))^(3/2)*cos(theta)*(2 + 3*gammat + ...
    6*n*gammat - exp(gammat + 2*n*gammat)*(2 + gammat + 2*n*gammat) +... 
   3*(1 + 2*n)^2*gammat*cos(theta))*sin(theta))/(2*(1 + 2*n)^2);
drhoA(2,1) = drhoA(1,2);
A = kron(rhoA.',speye(length(rhoA)))+kron(speye(length(rhoA)),rhoA);
drhot = drhoA.';
F = -full( 2*drhot(:).' * ( A \ drhoA(:) ) );