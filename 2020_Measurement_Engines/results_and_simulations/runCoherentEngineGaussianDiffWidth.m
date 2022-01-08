% function runCoherentEngineGaussian(nc)

% function getSteadyEngineXbasis(dimp,nh,nc,kh,kc,G,wp,ws,g)
%test measurement engine with x projectors, and HO in xbasis
clear;
Nx = 3;
wp = 1;
ws = 100;
nh = 1;
kh = 0.001;
kc = 0.1;
G = 0.1;
g = 2.5;
delta = g*g/wp;%[-5:0.5:-1 -0.8:0.2:0.8 0.825:0.025:1.3 1.4:0.2:3 3.5:0.5:5]*g*g/wp;
xeff = g/wp/sqrt(2);
mu = -xeff;
wd = ws-delta;
sigma = xeff*[0.05:0.05:1.5];
for nc = 1
    for i = 1:length(sigma)
        i
        %% Operators for qubit
        sz = [-1 0; 0 1];
        sp = [0 0; 1 0];
        sm = sp';
        sx = sp+sm;
        % sy = 1i*(sm-sp);
        idQ = eye(2);
        %% Operators for HO
        %hbar = 1
        %mass m=1
        %x = (a+adg)/sqrt(2), so unit should be sqrt(hbar/m/wp)
        %p = (a-adg)/1i/sqrt(2)
        %x interval from -xmax to xmax
        %choice critical! if xmax too high, then wrong results for thermal nbar at
        %low dimension dp. If too low, then probably wrong results for high nbar...
        dx = xeff/Nx;
        xmax = floor(10/dx)*(dx);
        x = -xmax:dx:xmax;
        dimp = length(x);
        p = linspace(-pi/dx,pi/dx,dimp); %fft ordering
        N = diag(0:dimp-1);
        a = diag(sqrt(1:dimp-1),1);
        adg = a';
        idHO = eye(dimp);
        %position operator in position repres.
        X = (diag(x.'));
        %FT basis trafo, from position to momentum repres.
        FTM = exp(-1i*p.'*x)/sqrt(dimp);
        % %check unitarity
        % fprintf( 'FT unitarity: %1.2e, %1.2e \n', norm(FTM'*FTM-eye(dimp)),norm(FTM*FTM'-eye(dimp)) );
        
        %momentum operator in position repres.
        P = FTM'*(diag(p.')*FTM);
        
        %shifted position operator due to interaction
        Xeff = kron(idQ,X) + g/wp/sqrt(2)*kron(sz,idHO);
        
        %driving function
        fx = exp(-(x-mu).^2/2/sigma(i)^2);
        fx = fx/max(fx);
        fx = diag(fx);
        
        %hamiltonian
        Hq = ws/2*kron(sz,idHO);
        Hp = wp/2*kron(eye(2),P*P) + wp/2*(Xeff*Xeff);
        Hfb = G*kron(sx,fx);
        Hrot = kron(sz,idHO);
        H = Hq + Hp + Hfb ;
        %coherent evolution
        spU = -1i* ( spLeftMultiply(H) - spRightMultiply(H) );
        spU2 = -1i* ( spLeftMultiply(Hrot) - spRightMultiply(Hrot));
        %hot bath dissipator
        spLh = kh*(nh+1)*(spLrMultiply(kron(sm,idHO)) - 0.5*spLeftMultiply(kron(sp*sm,idHO)) - 0.5*spRightMultiply(kron(sp*sm,idHO))) ...
            + kh*nh*(spLrMultiply(kron(sp,idHO)) - 0.5*spLeftMultiply(kron(sm*sp,idHO)) - 0.5*spRightMultiply(kron(sm*sp,idHO)));
        %cold bath dissipator
        A = (Xeff+1i*kron(idQ,P))/sqrt(2);
        Adg = A';
        spLc = kc*(nc+1)*(spLrMultiply(A) - 0.5*spLeftMultiply(Adg*A) - 0.5*spRightMultiply(Adg*A)) ...
            + kc*nc*(spLrMultiply(Adg) - 0.5*spLeftMultiply(A*Adg) - 0.5*spRightMultiply(A*Adg));
        %total
        spL = spLh+spLc+spU;
        spLtot = spL-spU2*wd/2;
        
        %% Solve for steady state
        sprhoSSAll = spnull(spLtot);
        sprhoSS = reshape(sprhoSSAll,dimp*2,dimp*2);
        % sprhoSS = sprhoSS.'; %need to transpose after reshape
        sprhoSS = sprhoSS/trace(sprhoSS);
        sprhog = sprhoSS(1:dimp,1:dimp);
        sprhoe = sprhoSS((1:dimp)+dimp,(1:dimp)+dimp);
        probg(i,:) = full(real(diag(sprhog)));
        probe(i,:) = full(real(diag(sprhoe)));
        
        %     rhoSSQ = ptrace(sprhoSS,2,dimp);
        plot(real(diag(sprhog))); hold on; plot(real(diag(sprhoe)));
        
        %% Get heat flows
        W(i) = 1i*G*wd*trace(kron(sm-sp,fx)*sprhoSS);
        Qh(i) = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*H);
        Qc(i) = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*H);
        %     W(i)
        % Qh1 = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hq);
        % Qc1 = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hq);
        % W1 = trace((G*diss(sprhoSS,kron(idQ,projP))+G*diss(sprhoSS,kron(sx,projM)))*Hq);
        % Qh2 = trace((kh*(nh+1)*diss(sprhoSS,kron(sm,idHO)) +kh*nh*diss(sprhoSS,kron(sp,idHO)))*Hp);
        % Qc2 = trace((kc*(nc+1)*diss(sprhoSS,A) +kc*nc*diss(sprhoSS,Adg))*Hp);
        % W2 = trace((G*diss(sprhoSS,kron(idQ,projP))+G*diss(sprhoSS,kron(sx,projM)))*Hp);
        % Qh3 = Qh1+Qh2;
        % Qc3 = Qc1+Qc2;
        % W3 = W1+W2;
        
        
    end
    clear sprhoSSAll sprhoSS sprhog sprhoe rhoSSQ spLtot spLc spLh spU spL spU2;
    filename = sprintf('Gaussian_delta1_nc%.3f.mat',nc);
    save(filename);
end