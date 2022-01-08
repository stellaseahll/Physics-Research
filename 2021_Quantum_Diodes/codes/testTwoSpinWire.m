%test script for twoSpinWire

%% run heat flows over coupling strength and spin numbers

Nvals = [3,3; 3,1; 1,3; 1,1];
sizeN = size(Nvals);
ks = 1e-3*[1,1];
Ts = [1,0.1];
ns = 1./ ( exp(1./Ts) - 1);
ws = [1,1];

gvals = logspace(-4,-1,101);


Q1 = zeros(sizeN(1),length(gvals));
Q2 = zeros(sizeN(1),length(gvals));

for k = 1:sizeN(1)
    

    tic
    s = twoSpinWire(Nvals(k,:),gvals(1),ks,ns,ws);
    s.steadyState();
    [q1,q2] = s.getHeatFlows();
    Q1(k,1) = q1;
    Q2(k,1) = q2;
    toc

    tic
    for j=2:length(gvals)

        s.g = gvals(j);
        s.initOperators();
        s.steadyState();
        [q1,q2] = s.getHeatFlows();
        Q1(k,j) = q1;
        Q2(k,j) = q2;

    end
    toc
end

figure();
semilogx(gvals,Q1,'LineWidth',1); xlabel('g'); ylabel('Qh'); 