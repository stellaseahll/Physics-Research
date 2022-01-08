clear;clc;
dimp = 30;
nc = 0.1;
G = 0.1;
wp = 1;
ws = 100;
nh = 1; %at ws
Th = ws/log(1/nh+1);
kh = 0.001;
kc = 0.1;
g = 2.5;
sigma = 1;

for idx1 = 1:length(g)
    delta = g(idx1)*g(idx1)/wp;
    wd = ws-delta;    
    C = zeros(dimp);
    for i = 0:dimp-1
        for j = i:dimp-1
            C(i+1,j+1) = getmDn(j,i,g(idx1)/wp);
            C(j+1,i+1) = C(i+1,j+1);
        end
    end
    sx = kron([0 1; 0 0], C) + kron([0 0; 1 0],C');
    sy = kron([0 1i; 0 0], C) -kron([0 0; 1i 0],C');
    sz = kron([-1 0; 0 1],eye(dimp));
    a = diag(sqrt(1:dimp-1),1);
    adg = a';
    b = kron(eye(2),a);
    bdg = kron(eye(2),adg);
    R = getR(dimp,g(idx1),wp,sigma*g(idx1)/sqrt(2)/wp);
    sxfx = kron([0 0;1 0],R) +  kron([0 1; 0 0],R');
    syfx = -kron([0 0; 1i 0],R) + kron([0 1i;0  0],R');
    Hq = ws/2*sz;
    Hp = wp*bdg*b;
    Hd = G*sxfx;
    Hrot = wd/2*sz;
    H = Hq + Hp + Hd - Hrot;
    
    spU = -1i* ( spLeftMultiply(H) - spRightMultiply(H) );
    X = cell(2*dimp+1);
    E = zeros(1,2*dimp+1);
    X{1} = kron([0 0; 1 0],diag(diag(C)));
    E(1) = ws;
    for k = 1:dimp
        E(1+k) = ws+wp*k;
        X{1+k} = sparse(kron([0 0; 1 0],diag(diag(C,-k),-k)));
        E(1+k+dimp) = ws-wp*k;
        X{1+k+dimp} = sparse(kron([0 0; 1 0],diag(diag(C,k),k)));
    end
    
    nh = 1./(exp(E./Th)-1);
    spLh = sparse(zeros(size(spU)));
    for i = 1:length(X)
        spLh = spLh + kh*(nh(i)+1)*(spLrMultiply(X{i}') - 0.5*spLeftMultiply(X{i}*X{i}') - 0.5*spRightMultiply(X{i}*X{i}'))...
            + kh*nh(i)*(spLrMultiply(X{i}) - 0.5*spLeftMultiply(X{i}'*X{i}) - 0.5*spRightMultiply(X{i}'*X{i}));
    end
    
    for idx2 = 1:length(nc)
        fprintf('(%d,%d)\n',idx1,idx2);
        spLc = kc*(nc(idx2)+1)*(spLrMultiply(b) - 0.5*spLeftMultiply(bdg*b) - 0.5*spRightMultiply(bdg*b)) ...
            + kc*nc(idx2)*(spLrMultiply(bdg) - 0.5*spLeftMultiply(b*bdg) - 0.5*spRightMultiply(b*bdg));
       
        spLtot = spLh+spLc+spU;
        %% Solve for steady state
        sprhoSSAll = spnull(spLtot);
        sprhoSS = reshape(sprhoSSAll,dimp*2,dimp*2);
        % sprhoSS = sprhoSS.'; %need to transpose after reshape
        sprhoSS = sprhoSS/trace(sprhoSS);
        sprhog = sprhoSS(1:dimp,1:dimp);
        sprhoe = sprhoSS((1:dimp)+dimp,(1:dimp)+dimp);
        rhoSSQ = ptrace(sprhoSS,2,dimp);
        pop{idx1,idx2} = real(diag(sprhoSS));
        
        W(idx1,idx2) = G*wd*trace(syfx*sprhoSS);
        % Qh = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*(H + Hrot));
        Qc(idx1,idx2) = trace((kc*(nc(idx2)+1)*diss(sprhoSS,b) +kc*nc(idx2)*diss(sprhoSS,bdg))*(H + Hrot));
        Qh(idx1,idx2) = 0;
        for i = 1:length(X)
            Qh(idx1,idx2) = Qh(idx1,idx2) + trace((kh*(nh(i)+1)*diss(sprhoSS,X{i}')+kh*(nh(i))*diss(sprhoSS,X{i}))...
                *(H + Hrot));
        end
    end
end

% save('Coherentwithfx_nc0.1.mat','nh','nc','kh','kc','g','G','W','Qh','Qc','ws','wp','pop');