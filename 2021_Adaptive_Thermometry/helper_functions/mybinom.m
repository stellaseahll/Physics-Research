function f = mybinom(nH,nTot,theta)

for i = 1:length(theta)
    f(i) = theta^nH * (1-theta)^(nTot-nH);
end
    