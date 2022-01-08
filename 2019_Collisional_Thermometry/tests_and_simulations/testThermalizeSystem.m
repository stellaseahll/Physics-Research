function rho = testThermalizeSystem(N) 
%thermalizing channel applied to system qubit 
%total state rho0 -> rho 
    d = 2^(N-1);
    E = rand();
    p = rand();
    tic;
    fprintf('create random state -- ');
    rho0 = rand(2^N);
    rho0 = rho0'*rho0;
    rho0 = rho0/trace(rho0);
    toc
    tic;
    fprintf('normal indexing -- ');
    %decay of coherences
    rho = sqrt(E)*rho0;
    %mixing of excited and ground populations
    rho(1:d,1:d) = (p+(1-p)*E) * rho0(1:d,1:d) + p*(1-E) * rho0(d+1:end,d+1:end);
    rho(d+1:end,d+1:end) = (1-p)*(1-E) * rho0(1:d,1:d) + (1-p+p*E) * rho0(d+1:end,d+1:end);
    toc
end
