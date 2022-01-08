%find max QC and g and chi for ohmic system at random temperatures

Ec = 1; %sets energy unit
%Ew = 1:0.5:10;
Ews = [1,2,3,4,5,6];
Ehs = Ews+Ec;

chiMin = 1e-5;
chiMax = 0.1;
gMin = 1e-4;
gMax = 1;

NT = 5000;
options = optimoptions('fmincon','Display','off','Algorithm','sqp');

for n=1:length(Ews)
    tic
    fprintf('Ew step %i -- ',n)
    Ew = Ews(n);
    Eh = Ehs(n);
        
    Tc = 1+rand(1,NT)*9;
    Th = Tc.*(1 + 0.9*rand(1,NT)*(Eh/Ec-1));
    Tw = Ew./(Eh./Th - Ec./Tc) .* (1 + 9*rand(1,NT));
    
    QC = zeros(1,NT);
    g = zeros(1,NT);
    chi = zeros(1,NT);
    keff = zeros(1,NT);
    
    for j=1:NT
        chiGuess = 0.1*6/( 2*Eh/(exp(Eh/Th(j))-1) + ...
            2*Ec/(exp(Ec/Tc(j))-1) + 2*Ew/(exp(Ew/Tw(j))-1) + Eh*2);
        if chiGuess>=chiMax
            chiGuess=chiMax/2;
        elseif chiGuess<=chiMin
            chiGuess = 2*chiMin;
        end
        [X,F] = fmincon(@(x)evalQ_modelOPX(x(1),x(2), [Th(j), Tc(j), Tw(j)], [Eh,Ec,Ew]), ...
            [0.1,chiGuess],[],[],[],[],[gMin,chiMin],[gMax,chiMax],[],options);
        QC(j) = -F;
        g(j) = X(1);
        chi(j) = X(2);
        keff(j) = g(j)*6/( 2*Eh/(exp(Eh/Th(j))-1) + ...
            2*Ec/(exp(Ec/Tc(j))-1) + 2*Ew/(exp(Ew/Tw(j))-1) + Eh*2);
    end
    fnam = sprintf('maxQCrandomT_Ec1_Ew%i_N5000_sqp.mat',Ew);
    save(fnam,'QC','g','chi','Tc','Tw','Th','Ec','Ew','Eh','keff','chiMax','chiMin','gMin','gMax');
    toc
end