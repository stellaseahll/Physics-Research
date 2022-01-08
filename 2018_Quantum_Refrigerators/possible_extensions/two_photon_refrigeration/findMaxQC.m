function D = findMaxQC(Ts,Es,optStep,gAcc)
g = linspace(0.001,1,optStep);
chi = linspace(0.001,0.1,optStep);
QHtmp = zeros(optStep);
QCtmp = zeros(optStep);
QWtmp = zeros(optStep);
gmax = g(end);
gmin = g(1);
chimin = chi(1);
chimax = chi(end);
while ((gmax-gmin)>gAcc)
    for j = 1:optStep
        for k = 1:optStep
            M = model3QubitsFridgeOld(g(k),Ts,chi(j)*Es,Es,'ohmic','x','p');
            M.findSS();
            QHtmp(j,k) = M.getTotalHeatFlow(1);
            QCtmp(j,k) = M.getTotalHeatFlow(2);
            QWtmp(j,k) = M.getTotalHeatFlow(3);
        end
    end
    [a,b] = find(QCtmp == max(max(QCtmp)));
    if (b==optStep)
       	gmin = max([gmin, g(b-1)]); 
        gmax = min([gmax, g(b)*2]);
    elseif (b==1)
        gmin = max([gmin, g(b)/2]);
        gmax = min([gmax, g(b+1)]);
    else
        gmin = g(b-1);
        gmax = g(b+1);
    end
    if (a==optStep)
       	chimin = max([chimin, chi(a-1)]); 
        chimax = min([chimax, chi(a)*2]);
    elseif (a==1)
        chimin = max([chimin, chi(a)/2]);
        chimax = min([chimax, chi(a+1)]);
    else
        chimin = chi(a-1);
        chimax = chi(a+1);
    end
    g = linspace(gmin,gmax,optStep);
%     gmin, gmax, chimin,chimax
    chi = linspace(chimin,chimax,optStep);
end

[a,b] = find(QCtmp == max(max(QCtmp)));
D = [chi(a), g(b), QCtmp(a,b), QCtmp(a,b)./QWtmp(a,b)];

