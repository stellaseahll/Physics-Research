function [Q,SS] = evalQ_modelOPX(g,chi,Ts,Es)

%internal functions
leftMultiply = @(A)(kron(eye(8),A)); %A*rho
rightMultiply = @(B)(kron(B.',eye(8))); %rho*B
lrMultiply = @(A)(kron(conj(A),A)); %A*rho*A'
lPLUSrMultiply = @(A)( kron(eye(8),A) + kron(A.',eye(8)) ); %A*rho + rho*A

%Hamilton
H = g*fliplr(eye(8)); % Sx x Sx x Sx
H = H + diag([0,Es(3),Es(2),Es(2)+Es(3),Es(1),Es(1)+Es(3),Es(1)+Es(2),Es(1)+Es(2)+Es(3)]);
[U , E] = eig(H,'vector');

X1 = [zeros(4),eye(4);eye(4),zeros(4)]; %S_x on qubit 1 (H)
X2 = blkdiag([0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0],[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]); %S_x on C
X3 = blkdiag([0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0],[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]); %S_x on W

%Transform into Eigenbasis for coefficients
expectX1 = U'*(X1*U);
expectX2 = U'*(X2*U);
expectX3 = U'*(X3*U);

%coarse graining time
deltat = 1;

%% bath 1
[r,c] = find(abs(expectX1)>10*eps);
Egaps = E(r) - E(c);
[Em,~,ic] = uniquetol(Egaps,10*eps); %careful, Em will be column vector here!
Egaps = Em(ic);
Nm = length(Em);

Emm = ones(Nm,1)*Em.';
W = (Emm.' + Emm)/2;
aW = abs(W);
omega = (Emm - Emm.')*deltat/2;
gamma = (aW.*(1./(exp(aW./Ts(1))-1) + heaviside(-W)));
gamma(aW<1e-3*Ts(1)) = Ts(1) + heaviside(-W(aW<1e-3*Ts(1))) ;
gamma(aW<10*eps) = Ts(1);
gamma = gamma.*exp(1i*omega).*sinc(omega/pi);

%no need to diagonalize, can do nondiag Lindblad double sum directly
[evec, ev] = eig(gamma,'vector');

K = cell(Nm,1);
for m=1:Nm
    K{m} = zeros(8);
    idx = find(Egaps==Em(m));
    for n=idx.'
        K{m} = K{m} + expectX1(r(n),c(n))*(U(:,r(n))*U(:,c(n))');
        %size(K{m})
    end
end

L1 = zeros(64);
for n=1:Nm
    LOpn = zeros(8);
    for k=1:Nm
        LOpn = LOpn + evec(k,n)*K{k};
    end
    L1 = L1 + chi*ev(n)*( lrMultiply(LOpn) - 0.5*lPLUSrMultiply(LOpn'*LOpn) );
%     L1 = L1 + chi*ev(n)*( lrMultiply(LOpn) ...
%         - 0.5*leftMultiply(LOpn'*LOpn) ...
%         - 0.5*rightMultiply(LOpn'*LOpn) );
end

%only nonzero elements of gamma required...
% [r,c] = find(abs(gamma)>10*eps);
% 
% L1 = zeros(64);
% for n=1:length(r)
%     KKn = K{c(n)}' * K{r(n)};
%     L1 = L1 + chi*gamma(r(n),c(n))*( leftMultiply(K{r(n)},8)*rightMultiply(K{c(n)}',8) ...
%         - 0.5*leftMultiply(KKn,8) ...
%         - 0.5*rightMultiply(KKn,8) );
% 
% end


%% bath 2
[r,c] = find(abs(expectX2)>10*eps);
Egaps = E(r) - E(c);
[Em,~,ic] = uniquetol(Egaps,10*eps); %careful, Em will be column vector here!
Egaps = Em(ic);
Nm = length(Em);

Emm = ones(Nm,1)*Em.';
W = (Emm.' + Emm)/2;
aW = abs(W);
omega = (Emm - Emm.')*deltat/2;
gamma = (aW.*(1./(exp(aW./Ts(2))-1) + heaviside(-W)));
gamma(aW<1e-3*Ts(2)) = Ts(2) + heaviside(-W(aW<1e-3*Ts(2))) ;
gamma(aW<10*eps) = Ts(2);
gamma = gamma.*exp(1i*omega).*sinc(omega/pi);

%no need to diagonalize, can do nondiag Lindblad double sum directly
[evec, ev] = eig(gamma,'vector');

K = cell(Nm,1);
for m=1:Nm
    K{m} = zeros(8);
    idx = find(Egaps==Em(m));
    for n=idx.'
        K{m} = K{m} + expectX2(r(n),c(n))*(U(:,r(n))*U(:,c(n))');
    end
end

L2 = zeros(64);
for n=1:Nm
    LOpn = zeros(8);
    for k=1:Nm
        LOpn = LOpn + evec(k,n)*K{k};
    end
    L2 = L2 + chi*ev(n)*( lrMultiply(LOpn) - 0.5*lPLUSrMultiply(LOpn'*LOpn) );
    %+ chi*ev(n)*( lrMultiply(LOpn) ...
    %    - 0.5*leftMultiply(LOpn'*LOpn) ...
    %    - 0.5*rightMultiply(LOpn'*LOpn) );
end

%only nonzero elements of gamma required...
% [r,c] = find(abs(gamma)>10*eps);
% 
% L2 = zeros(64);
% for n=1:length(r)
%     KKn = K{c(n)}' * K{r(n)};
%     L2 = L2 + chi*gamma(r(n),c(n))*( leftMultiply(K{r(n)},8)*rightMultiply(K{c(n)}',8) ...
%         - 0.5*leftMultiply(KKn,8) ...
%         - 0.5*rightMultiply(KKn,8) );
% end

%% bath 3
[r,c] = find(abs(expectX3)>10*eps);
Egaps = E(r) - E(c);
[Em,~,ic] = uniquetol(Egaps,10*eps); %careful, Em will be column vector here!
Egaps = Em(ic);
Nm = length(Em);

Emm = ones(Nm,1)*Em.';
W = (Emm.' + Emm)/2;
aW = abs(W);
omega = (Emm - Emm.')*deltat/2;
gamma = (aW.*(1./(exp(aW./Ts(3))-1) + heaviside(-W)));
gamma(aW<1e-3*Ts(3)) = Ts(3) + heaviside(-W(aW<1e-3*Ts(3))) ;
gamma(aW<10*eps) = Ts(3);
gamma = gamma.*exp(1i*omega).*sinc(omega/pi);

%no need to diagonalize, can do nondiag Lindblad double sum directly
[evec, ev] = eig(gamma,'vector');

K = cell(Nm,1);
for m=1:Nm
    K{m} = zeros(8);
    idx = find(Egaps==Em(m));
    for n=idx.'
        K{m} = K{m} + expectX3(r(n),c(n))*(U(:,r(n))*U(:,c(n))');
    end
end

L3 = zeros(64);
for n=1:Nm
    LOpn = zeros(8);
    for k=1:Nm
        LOpn = LOpn + evec(k,n)*K{k};
    end
    L3 = L3 + chi*ev(n)*( lrMultiply(LOpn) - 0.5*lPLUSrMultiply(LOpn'*LOpn) );
    %+ chi*ev(n)*( lrMultiply(LOpn) ...
    %    - 0.5*leftMultiply(LOpn'*LOpn) ...
    %    - 0.5*rightMultiply(LOpn'*LOpn) );
end

%only nonzero elements of gamma required...
% [r,c] = find(abs(gamma)>10*eps);
% L3 = zeros(64);
% for n=1:length(r)
%     KKn = K{c(n)}' * K{r(n)};
%     L3 = L3 + chi*gamma(r(n),c(n))*( leftMultiply(K{r(n)},8)*rightMultiply(K{c(n)}',8) ...
%         - 0.5*leftMultiply(KKn,8) ...
%         - 0.5*rightMultiply(KKn,8) );
% end

%% find SS and heat flow TO cold bath
SS = null(L1+L2+L3 -1i*leftMultiply(H) + 1i*rightMultiply(H));
SS = SS/sum(SS(1:9:64));
%size(SS)
% rSS = SS(:,1);
% SS = leftMultiply(H)*( L2*rSS );
% Q = -real( sum(SS(1:9:64)) );

%rSS = SS(:,1);
%Q = -real(trace(H * reshape( L2*SS,[8,8] )));
Q = -real(sum(conj(H(:)) .* (L2*SS) )); %same as trace if H hermitian
%rSS = reshape(rSS,[8,8]);



% for n=1:length(r)
%     Egap = E(r(n))-E(c(n));
%     idx = find(abs(Em-Egap)<eps*10);
%     if (isempty(idx))
%         Em = [Em Egap];
%         K{end+1} = (U(:,r(n))*U(:,c(n))')*expectX;
%     else
%         K{idx} = K{idx}+(obj.U(:,r(n))*obj.U(:,c(n))')*expectX;
%     end
% end

% 
% for j = 1:length(Em)
%     LOps{end+1} = zeros(8);
%     for k = 1:length(Em)
%         LOps{end} = LOps{end} + evec(k,j)*K{k};
% 
%     end
%     LOps{end} = sqrt(chi*ev(j))*LOps{end};
% end



