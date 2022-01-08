T = 1;
bathType = 'ohmic';
g=0.001;
w1 = 1;
w2 = 1;
Deltat = 1;
Omega = w1+w2;
Delta = w1-w2;
H = zeros(4);
H(1,1) = Omega;
H(2:3,2:3) = [Delta g; g -Delta];
H(4,4) = -Omega;
[v, E] = eig(H);
E = diag(E);
% raising and lowering op on 1st qubit
omega(1) = E(2)-E(4);
omega(2) = E(3)-E(4);
omega(3) = E(4)-E(2);
omega(4) = E(4)-E(3);
x = kron([0 1; 1 0],eye(2));
Lg = cell(4,1);
Lg{1} = (v(:,2)'*x*v(:,4))*v(:,2)*(v(:,4))' +  (v(:,1)'*x*v(:,3))*v(:,1)*(v(:,3))';
Lg{2} = (v(:,3)'*x*v(:,4))*v(:,3)*(v(:,4))' +  (v(:,1)'*x*v(:,2))*v(:,1)*(v(:,2))';
Lg{3} = Lg{1}';
Lg{4} = Lg{2}';
nbar = @(x) 1/(exp(x/T)-1);
switch bathType
    case 'ohmic'
        G = @(x) x;
    case 'flat'
        G = @(x) 1;
end
OmCoeff = @(x) exp(1i*x)*sinc(x);
% frequency matrices
for i = 1:4
    for j = 1:4
        nu = (omega(i) + omega(j))/2;
        Om = (omega(i) - omega(j));
        if (abs(nu)<10*eps)
            Gamma(i,j) = T*OmCoeff(Om*Deltat/2);
            continue;
        end
        Gamma(i,j) = G(abs(nu))*(nbar(abs(nu))+heaviside(-nu))*OmCoeff(Om*Deltat/2);
    end
end

%% Find Liouville coarsegrained
[evec, ev] = eig(Gamma,'vector');
Lp = cell(4,1);
Lvp = sparse(zeros(16));
for j = 1:4
    Lp{j} = sparse(zeros(4));
    for k = 1:4
        Lp{j} = Lp{j} + sparse(evec(k,j)*Lg{k});
    end
    Lvp = Lvp + ev(j)*(spLrMultiply(Lp{j})-0.5*spLeftMultiply(Lp{j}'*Lp{j})-0.5*spRightMultiply(Lp{j}'*Lp{j}));
end

%% Liouville global
Lvg = sparse(zeros(16));
for j = 1:4
    nu = omega(j);
    k =  G(abs(nu))*(nbar(abs(nu))+heaviside(-nu));
    Lvg = Lvg + k*(spLrMultiply(Lg{j})-0.5*spLeftMultiply(Lg{j}'*Lg{j})-0.5*spRightMultiply(Lg{j}'*Lg{j}));
end
nbar = @(x) 1/(exp(x/T)-1);
switch bathType
    case 'ohmic'
        G = @(x) x;
    case 'flat'
        G = @(x) 1;
end
OmCoeff = @(x) exp(1i*x)*sinc(x);
%% Liouville local
kup =  G(abs(w1))*(nbar(abs(w1)));
kdown =  G(abs(w1))*(nbar(abs(w1))+1);
down = sparse(kron([0 0; 1 0],eye(2)));
up = down';
Lvl = kup*(spLrMultiply(up)-0.5*spLeftMultiply(up'*down)-0.5*spRightMultiply(up'*down))+...
    kdown*(spLrMultiply(down)-0.5*spLeftMultiply(down'*up)-0.5*spRightMultiply(down'*up));

