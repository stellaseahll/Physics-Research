function J = getSpinOp(dim)

j = (dim-1)/2;
m = (-j:j).';
J{3} = diag(m); %Jz
t = sqrt((j-m).*(j+m+1));
J{4} = diag(t(1:end-1),-1); %Jp
J{5} = J{4}'; %Jm
J{1} = 0.5*(J{4} + J{5}); %Jx
J{2} = 0.5*1i*(J{4} - J{5}); %Jy
k = zeros(dim,1);
k(m<0)=1;
J{6} = diag(k); %Proj-
k = zeros(dim,1);
k(m>=0)=1;
J{7} = diag(k); %Proj+