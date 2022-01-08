function K = bestFit(x,alpha,Wdata)

K = sum((Wdata - (x(1)+x(2)*log(alpha))).^2);