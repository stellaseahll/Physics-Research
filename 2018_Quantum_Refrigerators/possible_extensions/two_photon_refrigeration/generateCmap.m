function map = generateCmap(x)
k = [min(min(x)); 0;max(max(x))]-min(min(x));
k = k/max(k);
map = [1 0 0; 1 1 1; 0 0 1];
map = interp1(k,map,linspace(0,1,255));