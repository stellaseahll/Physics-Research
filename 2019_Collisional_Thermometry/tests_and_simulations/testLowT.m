%What's the best 1,2-ancilla thermometry strategy at low T?
%What to expect:
% -coupling g shouldn't be too small. But is full swap always optimal?
% -thermalisation time shouldn't be too large, otherwise thFI...

%% brute force approach: random parameters and compare QFI to thFI
% Trend: ground-state ancillas & full swaps are always best, 
% but the improvement over thermal FI decreases with T and 
% shifts to greater optimal gammat. 
% That is, in the limit of low T, ful thermalisation and th. FI always wins.

N = 1e5;
T = 0.3;
nbar = 1/(exp(1/T)-1);
dnbar = 1e-4; %relative difference

p1th = nbar/(2*nbar+1);
Fth = 1/(nbar*(nbar+1)*(2*nbar+1)^2);
%Fth2 = (p1th-p1th^2)/( (exp(1/T)/(exp(1/T)-1)^2)^2 );

gt = rand(1,N)*pi;
gammat = rand(1,N)*10;
pA = rand(1,N);
rA = rand(1,N).*sqrt(pA.*(1-pA)); %coherence, phase shouldn't play a role

F1 = zeros(1,N);

tic 
parfor i = 1:N
    rhoA = [pA(i), rA(i); rA(i), 1-pA(i)];
    s = spSinglePassXY(gt(i),gammat(i),nbar,dnbar*nbar,rhoA,1,1);
    s.alg = 3;
    F1(i) = s.getAllFish();
end
toc

[maxF,nmaxF] = max(F1/Fth);
fprintf('max F1/Fth = %1.4f at gt/pi=%1.2f, gammat=%1.2f, pA=%1.2f, rA=%1.2f \n', ...
    maxF, gt(nmaxF)/pi, gammat(nmaxF), pA(nmaxF), rA(nmaxF));

%% OK try 2-ancilla brute force, maybe better due to correlations...
% Doesn't seem to improve anything. At higher T > 0.1 it is still always
% best to use GS probes and full swap to beat SQL. At lower T, SQL is not
% reached, and GS with full swap still achieves max. Nothing new...

N = 1e5;
T = 0.5;
nbar = 1/(exp(1/T)-1);
dnbar = 1e-4; %relative difference

p1th = nbar/(2*nbar+1);
Fth = 1/(nbar*(nbar+1)*(2*nbar+1)^2);
%Fth2 = (p1th-p1th^2)/( (exp(1/T)/(exp(1/T)-1)^2)^2 );

gt = rand(1,N)*pi;
gammat = rand(1,N)*10;
pA = rand(1,N);
rA = rand(1,N).*sqrt(pA.*(1-pA)); %coherence, phase shouldn't play a role

F1 = zeros(1,N);
F2 = zeros(1,N);

tic 
parfor i = 1:N
    rhoA = [pA(i), rA(i); rA(i), 1-pA(i)];
    s = spSinglePassXY(gt(i),gammat(i),nbar,dnbar*nbar,rhoA,2,1);
    s.alg = 3;
    F = s.getAllFish();
    F1(i) = F(1);
    F2(i) = F(2);
end
toc

[maxF,nmaxF] = max(F1/Fth);
fprintf('max F1/Fth = %1.4f at gt/pi=%1.2f, gammat=%1.2f, pA=%1.2f, rA=%1.2f \n', ...
    maxF, gt(nmaxF)/pi, gammat(nmaxF), pA(nmaxF), rA(nmaxF));

[maxF,nmaxF] = max(F2/2/Fth);
fprintf('max F2/2Fth = %1.4f at gt/pi=%1.2f, gammat=%1.2f, pA=%1.2f, rA=%1.2f \n', ...
    maxF, gt(nmaxF)/pi, gammat(nmaxF), pA(nmaxF), rA(nmaxF));

%% OK now weak coupling 2-ancilla brute force. Maybe here something interesting...
% If we allow for strong coupling, full swap always performs best using GS
% probes and no correlations. Boring! So let's see whether there is more
% interesting behavior for 1,2 probes if we fix weak coupling!
%
% Once again, nothing interesting. At low T, there is almost no 2-probe
% enhancement, and max QFI is achieved with GS probes. Only at higher T, we
% find enhancements, and max QFI is obtained with pure states close to GS,
% but with a bit of coherence. 
% This makes sense somehow: At low T, the relaxation rate to the bath is
% almost T-independent, so by looking at correlations between 2 subsequent
% probes, we don't learn much about T from the decoherence in between. 
% --> Low T is simply boring...

N = 1e5;
T = 2.5;
nbar = 1/(exp(1/T)-1);
dnbar = 1e-4; %relative difference

p1th = nbar/(2*nbar+1);
Fth = 1/(nbar*(nbar+1)*(2*nbar+1)^2);
%Fth2 = (p1th-p1th^2)/( (exp(1/T)/(exp(1/T)-1)^2)^2 );

gt = 0.1*pi; %fixed weak value
gammat = rand(1,N)*10;
pA = rand(1,N);
rA = rand(1,N).*sqrt(pA.*(1-pA)); %coherence, phase shouldn't play a role

F1 = zeros(1,N);
F2 = zeros(1,N);

tic 
parfor i = 1:N
    rhoA = [pA(i), rA(i); rA(i), 1-pA(i)];
    s = spSinglePassXY(gt,gammat(i),nbar,dnbar*nbar,rhoA,2,1);
    s.alg = 3;
    F = s.getAllFish();
    F1(i) = F(1);
    F2(i) = F(2);
end
toc

[maxF1,nmaxF] = max(F1/Fth);
fprintf('max F1/Fth = %1.4f at gt/pi=%1.2f, gammat=%1.2f, pA=%1.2f, rA=%1.2f \n', ...
    maxF1, gt/pi, gammat(nmaxF), pA(nmaxF), rA(nmaxF));

[maxF2,nmaxF] = max(F2/2/Fth);
fprintf('max F2/2Fth = %1.4f at gt/pi=%1.2f, gammat=%1.2f, pA=%1.2f, rA=%1.2f \n', ...
    maxF2, gt/pi, gammat(nmaxF), pA(nmaxF), rA(nmaxF));

figure(); plot(F1/Fth,F2/2/Fth,'.',[0 maxF1], [0 maxF1],'k--')