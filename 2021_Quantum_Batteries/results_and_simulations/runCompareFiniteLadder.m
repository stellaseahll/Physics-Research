%fix system as qubit

clear;clc;
invT = -logspace(-1,2,50);
dB = 2:100;
for j = 1:length(invT)
    j
    for k = 1:length(dB)
        fprintf('(%d,%d)\n',j,k)
        s = chargeBattery(invT(j),2,dB(k),2);
        s.runSameInteraction(100);
        Etype2{j,k} = s.getEnergy();
        E2type2{j,k} = s.getE2();
        p = chargeBattery(invT(j),2,dB(k),3);
        p.runSameInteraction(100);
        Etype3{j,k} = p.getEnergy();
        E2type3{j,k} = p.getE2();
        tmp2 = diff(Etype2{j,k});
        if abs((tmp2(2)-tmp2(1))/tmp2(2)+tmp2(1))>1e-12 %to check if boundary is reached.
            dEtype2(j,k) = tmp2(1);
        end
        tmp3 = diff(Etype3{j,k});
        if abs((tmp3(2)-tmp3(1))/tmp3(2)+tmp3(1))>1e-12 %to check if boundary is reached.
            dEtype3(j,k) = tmp3(1);
        end
    end
end

