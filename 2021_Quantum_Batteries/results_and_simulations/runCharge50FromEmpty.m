%Compute how many steps/time needed to achieve X% battery charge, as
%function of q,theta for coherent charging, can compare to "ideal"
%classical strategy of using q=0 and theta=pi/2 (deterministic charging)
N=200;
Ngoal = 0.1*N;

q = linspace(0,0.5,101);
theta = (1:50)/50*pi/2;
K50 = zeros(length(q),length(theta));

for j = 1:length(q)
    fprintf('row %i/%i -- ',j,length(q));
    tic
    for k = 1:length(theta)
        K50(j,k) = chargeStepsFromEmpty(q(j),theta(k),N,Ngoal);
    end
    toc
end

T50 = K50.*theta;