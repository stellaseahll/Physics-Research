function [x f] = myMinSearch(fun,a,b,tol)
%not global minimum
%one minimum between [a,b], otherwise not guaranteed to be global
xstart =a; 
xend = b;
fstart = fun(xstart);
fend = fun(xend);
dx = b-a;
while (dx>tol)
    xmid = (xstart+xend)/2;
    if (fstart<=fend)
        xend = xmid;
        fend = fun(xend);
    else
        xstart = xmid;
        fstart = fun(xstart);
    end
    dx = xend-xstart;
end

if fstart<fend
    f = fstart;
    x = xstart;
else
    f = fend;
    x = xend;
end