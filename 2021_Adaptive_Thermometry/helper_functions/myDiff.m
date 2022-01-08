function dydx = myDiff(dx,y)
dydx = [(y(2)-y(1))/dx (y(3:end)-y(1:end-2))/2/dx (y(end)-y(end-1))/dx];