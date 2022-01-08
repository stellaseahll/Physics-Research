function rho = prepareBell(k)
rho = zeros(4);
switch k
    case 1 %gg+ee
        rho(1,1)=0.5;
        rho(1,4)=0.5;
        rho(4,1)=0.5;
        rho(4,4)=0.5;
    case 2 %gg-ee
        rho(1,1)=0.5;
        rho(1,4)=-0.5;
        rho(4,1)=-0.5;
        rho(4,4)=0.5;
    case 3 %ge+eg
        rho(2,2)=0.5;
        rho(3,3)=0.5;
        rho(2,3)=0.5;
        rho(3,2)=0.5;
    case 4 %ge-eg
        rho(2,2)=0.5;
        rho(3,3)=0.5;
        rho(2,3)=-0.5;
        rho(3,2)=-0.5;
end