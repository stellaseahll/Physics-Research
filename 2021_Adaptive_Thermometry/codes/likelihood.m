function f = likelihood(nH,nTot,theta)

for i = 1:length(theta)
    f(i) = theta(i)^nH * (1-theta(i))^(nTot-nH);
end
    