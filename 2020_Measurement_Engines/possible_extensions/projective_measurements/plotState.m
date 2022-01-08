function plotState(rho,mode)
    if (nargin == 1)
        mode = 1;
    end
    mz = (rho(1) - rho(4));
    mx = real(rho(1,2))*2;
    my = imag(rho(1,2))*2;
    if (mode == 2)
        plot3([0,mx],[0,my],[0,mz],'linewidth',2);
    else 
        plot([0,mx],[0,mz],'linewidth',2);
    end
end
