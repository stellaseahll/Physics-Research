function myhist(y,x)
    y = sort(y)
    miny = y[1];
    maxy = y[end];
    dy = (maxy-miny)/(x*1.0)
    edges = zeros(x)
    xo = zeros(x)
    for i = 1:x
        edges[i] = miny + i*dy;
        xo[i] = miny - 0.5*dy + i*dy;
    end
    # edges = linspace(miny,maxy,x+1);
    idx = 1;
    pdf = zeros(x+1)
    pdf[end] = length(y)
    for k = 1:length(y)
        while (y[k]>edges[idx] && idx<x)
            pdf[idx+1]= k-1
            idx += 1
        end
    end
    no = (pdf[2:end]-pdf[1:end-1])/length(y)
    return (no,xo)
end