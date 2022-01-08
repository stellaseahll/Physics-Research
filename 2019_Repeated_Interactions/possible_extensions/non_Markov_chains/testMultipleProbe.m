% 2018-09-10
% Test for multiple probe repeated interactions and memory effects
% Set-up: n probes interacting with system for time tau, trace out 1st bath
% qubit and replace with another bath qubit.
clear;clc;
%% Parameters
omegaS = 1.1;
omegaP = 1.1;
T = 1;
g = 0.1;
% tau = 10; %interaction time
% dt = 0; %free evolution time
N = 4; %number of probes
rhoP = [1 exp(-omegaP/T)];
rhoP = diag(rhoP)/sum(rhoP);
rhoS = [1 0; 0 0];
%% Operators
sz = [-1 0; 0 1];
sx = [0 1; 1 0];
sp = [0 0; 1 0];
sm = sp';
swap = zeros(4);
swap(1,1) = 1;
swap(2,3) = 1;
swap(3,2) = 1;
swap(4,4) = 1;
dt = 0.1;
tstep = 0:dt:2500; 
tau = [0.1 1 10 100]; % Position it traces over
for n = 1:N
    n
    for j = 1:length(tau)
    S = sparse(kron(swap,eye(2^(n-1))));
    H0 = kron(sz*omegaS,eye(2^n));
    Hint = zeros(2^(n+1));
    rhoNP = 1;
    for i = 1:n
        A = eye(2^i);
        B = eye(2^(n-(i+1)+1));
        C = eye(2^(i-1));
        H0 = H0 + omegaP*kron(A,kron(sz,B));
        Hint = Hint + g*kron(sp,kron(C,kron(sm,B))) + g*kron(sm,kron(C,kron(sp,B)));
        rhoNP = kron(rhoNP,rhoP);
    end
    rho = sparse(kron(rhoS,rhoNP));
    rho2 = rho;
    H = sparse(H0 + Hint);
%     U0 = expm(-1i*H0*dt);
    U = expm(-1i*H*dt);
%     Utotal = S*U*U0;
    sz1 = sparse(kron(sz,eye(2^n)));
%     sz2 = sparse(kron(eye(2),kron(sz,eye(2^(n-1)))));
        for i = 1:length(tstep)
            ez{j,n}(i) = sum(sum(sz1.*rho));
            ez2{j,n}(i) = sum(sum(sz1.*rho2));
            rho = U*rho*U';
            rho2 = U*rho2*U';
            if (mod(tstep(i),tau(j))==0)
                [~,rho] = ptrace(S*rho*S',2,2^n);
                rho = sparse(kron(rho,rhoP));
                rho2 = ptrace(rho2,2,2^n);
                rho2 = sparse(kron(rho2,rhoNP));
            end
        end
        save('data.mat');
    end
    %     plot(1:tstep,real(ez),[1 tstep],[1 1]*(rhoP(end)-rhoP(1))); hold on;
end

% omegaS = 1.1;
% omegaP = 1.1;
% T = 1;
% g = 0.01;
% % tau = 10; %interaction time
% % dt = 0; %free evolution time
% N = 6; %number of probes
% rhoP = [1 exp(-omegaP/T)];
% rhoP = diag(rhoP)/sum(rhoP);
% rhoS = [1 0; 0 0];
% %% Operators
% sz = [-1 0; 0 1];
% sx = [0 1; 1 0];
% sp = [0 0; 1 0];
% sm = sp';
% swap = zeros(4);
% swap(1,1) = 1;
% swap(2,3) = 1;
% swap(3,2) = 1;
% swap(4,4) = 1;
% dt = 0.1;
% tstep = 0:dt:40000; 
% tau = [0.1 0.2 0.5 1 2 5 10 20 50 100 200 500]; % Position it traces over
% for n = 1:N
%     n
%     for j = 1:length(tau)
%     S = sparse(kron(swap,eye(2^(n-1))));
%     H0 = kron(sz*omegaS,eye(2^n));
%     Hint = zeros(2^(n+1));
%     rhoNP = 1;
%     for i = 1:n
%         A = eye(2^i);
%         B = eye(2^(n-(i+1)+1));
%         C = eye(2^(i-1));
%         H0 = H0 + omegaP*kron(A,kron(sz,B));
%         Hint = Hint + g*kron(sp,kron(C,kron(sm,B))) + g*kron(sm,kron(C,kron(sp,B)));
%         rhoNP = kron(rhoNP,rhoP);
%     end
%     rho = sparse(kron(rhoS,rhoNP));
%     rho2 = rho;
%     H = sparse(H0 + Hint);
% %     U0 = expm(-1i*H0*dt);
%     U = expm(-1i*H*dt);
% %     Utotal = S*U*U0;
%     sz1 = sparse(kron(sz,eye(2^n)));
% %     sz2 = sparse(kron(eye(2),kron(sz,eye(2^(n-1)))));
%         for i = 1:length(tstep)
%             ez{j,n}(i) = sum(sum(sz1.*rho));
%             ez2{j,n}(i) = sum(sum(sz1.*rho2));
%             rho = U*rho*U';
%             rho2 = U*rho2*U';
%             if (mod(tstep(i),tau(j))==0)
%                 [~,rho] = ptrace(S*rho*S',2,2^n);
%                 rho = sparse(kron(rho,rhoP));
%                 rho2 = ptrace(rho2,2,2^n);
%                 rho2 = sparse(kron(rho2,rhoNP));
%             end
%         end
%         save('data2.mat');
%     end
% 
%     %     plot(1:tstep,real(ez),[1 tstep],[1 1]*(rhoP(end)-rhoP(1))); hold on;
% end

