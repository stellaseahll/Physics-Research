function [x f] = myMinSearch(fun,a,b,tol)
%not global minimum
%one minimum between [a,b], otherwise not guaranteed to be global
xstart =a; 
xend = b;
dx = b-a;
while (dx>tol)
    xmid = (xstart+xend)/2;
    if (fun(xstart)<=fun(xend))
        xend = xmid;
    else
        xstart = xmid;
    end
    dx = xend-xstart;
end

if fun(xstart)<fun(xend)
    f = fun(xstart);
    x = xstart;
else
    f = fun(xend);
    x = xend;
end