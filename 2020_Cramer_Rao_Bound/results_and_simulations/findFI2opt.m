function [F rhoA]= findFI2opt(theta,gammat,n)

% rho = (cos(theta/2) ; sin(theta/2))*(cos(theta/2) ; sin(theta/2))';
% s = spsinglePassXYZ(pi/2,0,gammat,n,1e-7,rho,2,1); s.alg = 1;
% f = s.getAllFish();
% F = -real(f(2));

L = 1-exp(-gammat*(2*n+1));

rhoA(1,1) = (1 + 2*n - L + (1 + 2*n)*(-1 + L)*cos(theta))^2/(4*(1 + 2*n)^2);
rhoA(1,2) = (sqrt(1 - L)*(1 + 2*n - L + (1 + 2*n)*(-1 + L)*cos(theta))*sin(2*theta))/(8 + 16*n);
rhoA(1,3) = (sqrt(1 - L)*((1 + 2*n)*(-1 + L) + (1 + 2*n - L)*cos(theta))*(-L + (1 + 2*n)*(-1 + L)*cos(theta))*sin(theta))/(4*(1 + 2*n)^2);
rhoA(1,4) =-(((-1 + L)*cos(theta)*(-L + (1 + 2*n)*(-1 + L)*cos(theta))*sin(theta)^2)/(4 + 8*n));
rhoA(2,2) = ((1 + 2*n)^2 - L^2 + (1 + 2*n)*(-1 + L)*cos(theta)*(2*L - (1 + 2*n)*(-1 + L)*cos(theta)))/(4*(1 + 2*n)^2);
rhoA(2,3) =((-1 + L)*cos(theta)*(-L + (1 + 2*n)*(-1 + L)*cos(theta))*sin(theta)^2)/(4 + 8*n);
rhoA(2,4) =(sqrt(1 - L)*(-L + (1 + 2*n)*(-1 + L)*cos(theta))*(-(1 + 2*n)*(-1 + L) + (1 + 2*n + L)*cos(theta))*sin(theta))/(4*(1 + 2*n)^2);
rhoA(3,3) = ((1 + 2*n)^2 - L^2 + (1 + 2*n)*(-1 + L)*cos(theta)*(2*L - (1 + 2*n)*(-1 + L)*cos(theta)))/(4*(1 + 2*n)^2);
rhoA(3,4) = (sqrt(1 - L)*(-1 - 2*n - L + (1 + 2*n)*(-1 + L)*cos(theta))*sin(2*theta))/(8 + 16*n);
rhoA(4,4) = (1 + 2*n + L - (1 + 2*n)*(-1 + L)*cos(theta))^2/(4*(1 + 2*n)^2);
rhoA(2,1) = rhoA(1,2);
rhoA(3,1) = rhoA(1,3);
rhoA(4,1) = rhoA(1,4);
rhoA(3,2) = rhoA(2,3);
rhoA(4,2) = rhoA(2,4);
rhoA(4,3) = rhoA(3,4);

drhoA(1,1)=-((exp(-2*(1 + 2*n)*gammat)*(-1 - 2*exp(gammat + 2*n*gammat)*n + cos(theta) + ...
    2*n*cos(theta))*(-1 + exp(gammat + 2*n*gammat) - gammat - 2*n*gammat + (1 + 2*n)^2*gammat*cos(theta)))/(1 + 2*n)^3);
drhoA(1,2)=((exp(-(1 + 2*n)*gammat))^(3/2)*(-2 - 3*gammat - 6*n*gammat - ...
   2*exp(gammat + 2*n*gammat)*(-1 + n*(1 + 2*n)*gammat) + 3*(1 + 2*n)^2*gammat*cos(theta))*sin(2*theta))/(8*(1 + 2*n)^2);
drhoA(1,3)=(1/(4*(1 + 2*n)^3))*exp(-2*(gammat + 2*n*gammat))*sqrt(exp(-(1 + 2*n)*gammat))*((1 + 2*n)*(2 + 5*(1 + 2*n)*gammat - ...
      exp(gammat + 2*n*gammat)*(2 + (3 + 6*n)*gammat)) + cos(theta)*(-4 - 10*(1 + 2*n*(2 + n*(3 + 2*n)))*gammat + 2*exp(2*(1 + 2*n)*gammat)*(-1 + ...
         n*(2 + gammat + 2*n*gammat)) + exp(gammat + 2*n*gammat)*(3*(2 + gammat) - 4*n*(1 + 3*n*gammat)) + (1 + 2*n)*(2 + 5*gammat + ...
         10*n*gammat + 2*exp(gammat + 2*n*gammat)*(-1 + 3*n*(1 + 2*n)*gammat))*cos(theta)))*sin(theta);
drhoA(1,4)=(exp(-2*(1 + 2*n)*gammat)*cos(theta)*(-1 - 2*gammat - 4*n*gammat + exp(gammat + 2*n*gammat)*(1 + gammat + 2*n*gammat) + ...
   2*(1 + 2*n)^2*gammat*cos(theta))*sin(theta)^2)/(2*(1 + 2*n)^2);
drhoA(2,2)=(exp(-2*(1 + 2*n)*gammat)*(-1 + exp(gammat + 2*n*gammat) + cos(theta) + 2*n*cos(theta))*(-1 + exp(gammat + 2*n*gammat) - ...
    gammat - 2*n*gammat + (1 + 2*n)^2*gammat*cos(theta)))/(1 + 2*n)^3;
drhoA(2,3)=(exp(-2*(1 + 2*n)*gammat)*cos(theta)*(1 + 2*gammat + 4*n*gammat - ...
   exp(gammat + 2*n*gammat)*(1 + gammat + 2*n*gammat) - 2*(1 + 2*n)^2*gammat*cos(theta))*sin(theta)^2)/(2*(1 + 2*n)^2);
drhoA(2,4)=(1/(4*(1 + 2*n)^3))*exp(-2*(1 + 2*n)*gammat)*sqrt(exp(-(1 + 2*n)*gammat))*((1 + 2*n)*(-2 - 5*(1 + 2*n)*gammat + ...
      exp(gammat + 2*n*gammat)*(2 + (3 + 6*n)*gammat)) + cos(theta)*(4 + 10*(1 + 2*n*(2 + n*(3 + 2*n)))*gammat + 2*exp(2*(1 + 2*n)*gammat)*(3 + gammat + ... 
         n*(2 + 3*gammat + 2*n*gammat)) - exp(gammat + 2*n*gammat)*(10 + 9*gammat + 4*n*(1 + 3*(2 + n)*gammat)) + (1 + 2*n)*(-2 - 5*(1 + 2*n)*gammat + ...
         2*exp(gammat + 2*n*gammat)*(1 + 3*(1 + n)*(1 + 2*n)*gammat))*cos(theta)))*sin(theta);
drhoA(3,3)=(exp(-2*(1 + 2*n)*gammat)*(-1 + exp(gammat + 2*n*gammat) + cos(theta) + 2*n*cos(theta))*(-1 + exp(gammat + 2*n*gammat) - gammat - ...
   2*n*gammat + (1 + 2*n)^2*gammat*cos(theta)))/(1 + 2*n)^3;
drhoA(3,4)=((exp(-(1 + 2*n)*gammat))^(3/2)*(-2 - 3*gammat - 6*n*gammat + 2*exp(gammat + ...
     2*n*gammat)*(1 + (1 + n)*(1 + 2*n)*gammat) + 3*(1 + 2*n)^2*gammat*cos(theta))*sin(2*theta))/(8*(1 + 2*n)^2);
drhoA(4,4)=-((exp(-2*(1 + 2*n)*gammat)*(-1 + 2*exp(gammat + 2*n*gammat)*(1 + n) + cos(theta) + ...
    2*n*cos(theta))*(-1 + exp(gammat + 2*n*gammat) - gammat - 2*n*gammat + (1 + 2*n)^2*gammat*cos(theta)))/(1 + 2*n)^3);
drhoA(2,1) = drhoA(1,2);
drhoA(3,1) = drhoA(1,3);
drhoA(4,1) = drhoA(1,4);
drhoA(3,2) = drhoA(2,3);
drhoA(4,2) = drhoA(2,4);
drhoA(4,3) = drhoA(3,4);
A = kron(rhoA.',speye(length(rhoA)))+kron(speye(length(rhoA)),rhoA);
drhot = drhoA.';
F = -full( 2*drhot(:).' * ( A \ drhoA(:) ) );
       
end