
minE = 0:100:5000;
Tmin = 1;
Tmax = 10;
NT = 300;
N = 1:450;
clear optE;
clear minfE;
for i = 1:length(N)
    i
    for j = 1:length(minE)-1
        [E(j),fE(j)] =  fmincon(@(x) findOptBoundLargeN(x,N(i),prior,T),(minE(j+1)-minE(j))/2,[],[],[],[],minE(j),minE(j+1),[],options);
    end
    optE(i) = E(fE==min(fE));
    minfE(i) = fE(fE==min(fE));
end

const= 1.00466430;
newBound2 = 1./(F+(1:1e4)*const);
