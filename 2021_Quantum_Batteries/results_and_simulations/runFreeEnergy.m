load('fig1NEW.mat');
for i = 1:length(q)
    sEff2(i,2:tend+1) = sEff(i,2:tend+1)./(1:tend);
    tEff2(i,2:tend+1) = tEff(i,2:tend+1)./(1:tend);
    sEEff2(i,2:tend+1) = sEEff(i,2:tend+1)./(1:tend);
    tEEff2(i,2:tend+1) = tEEff(i,2:tend+1)./(1:tend);
end
theta = pi/2;
T = logspace(-2,1,100);
for i = 1:length(q)
    i
    sEff = sEff2(i,251);
    tEff = tEff2(i,251);
    for j = 1:length(T)
        Erg(i,j) = 1-2*q(i);
        qth = exp(0.5/T(j))/(exp(-0.5/T(j))+exp(0.5/T(j)));
        if ((1-qth)<eps)
            Fth(i,j) = (1-qth) + T(j)*(qth*log(qth));
        else
            Fth(i,j) = (1-qth) + T(j)*(qth*log(qth) + (1-qth)*log(1-qth));
        end
        ErgQ(i,j) = 1-q(i);
        if q(i)<eps
            F(i,j) = (1-q(i)) + T(j)*((1-q(i))*log(1-q(i)));
            FQ(i,j) = (1-q(i));
        else
            F(i,j) = (1-q(i)) + T(j)*(q(i)*log(q(i))+ (1-q(i))*log(1-q(i)));
            FQ(i,j) = (1-q(i));
        end
        sEff_k250(i,j) = sEff2(i,251)*Erg(i,j)/(F(i,j)-Fth(i,j))*250;
        tEff_k250(i,j) = tEff2(i,251)*ErgQ(i,j)/(FQ(i,j)-Fth(i,j))*250;
    end
end
