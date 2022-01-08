Tw = 8;
Tc = 1;

Nh = 51;
Th = Tc + (1:Nh)/Nh*(Tw-Tc);

Es = [5,1,4];

gs = zeros(1,Nh);
chis = zeros(1,Nh);
QCs = zeros(1,Nh);
etas = zeros(1,Nh);


for j=1:Nh
    tic
    %Th(j)
    fprintf('step %i -- ',j)
    D = findMaxQC([Th(j),Tc,Tw],Es,5,1e-4);
    chis(j) = D(1);
    gs(j) = D(2);
    QCs(j) = D(3);
    etas(j) = D(4);
    toc
end