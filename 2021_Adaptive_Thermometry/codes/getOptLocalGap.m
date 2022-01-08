function E = getOptLocalGap(Np)

global N
N = 2^Np;
E = fsolve(@(x) findGap(x),10);

end


function y = findGap(x)
global N

y = (N-1)*(x+2)/(x-2)-exp(x);
end