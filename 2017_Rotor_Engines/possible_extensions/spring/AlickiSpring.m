% Test script Alicki's spring

% x shall be coordinate of spring, in units of valve opening length L.
% x=0 shall be spring equilibrium IF no push by working mode (i.e. at T=0)
% valves such that radiation pressure pushes towards x>0
Nx = 201;
x = linspace(-3,3,Nx);

%valve coupling functions (corresponds to f_H^2 and f_C^2 in our engine, btw 0 and 1)
% fH = @(x)( 1 + (abs(x-1)-1).*heaviside(x).*heaviside(2-x) );
% fC = @(x)( (abs(x)-1).*heaviside(1-abs(x)) + 1 );

%any reason why we dont set it to be symmetric about x = 0? Since x=0 is
%the eqm position
fH = @(x)( 1 + (abs(x-0.5)-1).*heaviside(x+0.5).*heaviside(1.5-x) );
fC = @(x)( (abs(x+0.5)-1).*heaviside(1-abs(x+0.5)) + 1 );



%Assuming linear potential, i.e. homogeneous force due to radiation
%pressure, i.e. V(x) = -hbar*g*n*x/L with n the work mode occupation.
% Ideal fast thermalization, spring sees instantaneous rad. press. force:
%  F(x) ~ hbar*g/L*nth(x), nth(x) = ( nH*fH(x) + nC*fC(x) )/( fH(x) + fC(x) )
% that is on top of linear restoration force of spring!
nH = 1;
nC = 0;
w = 1; %spring eigenfrequency (if 1, then time in units 1/w)
gamma = 0.01; %possible damping of spring

%eff mass translation rate r=sqrt(hbar*g/(m*L^2)), determines spring stiffness VS pressure
r = 0.5;
%if r << w then free oscillation regime with small gain per cycle, can
%integrate over instantaneous rad press force to get per cycle gain, etc.

%instantaneous rad. press. accel. if fast thermalization
aRP = @(x)( r^2*( nH*fH(x) + nC*fC(x) )./( fH(x) + fC(x) ) );

%% test plots

%test plot coupling functions
figure(); 
plot(x,fC(x),x,fH(x));
grid on
title('coupling functions fC (blue), fH (red)')

%test plot accelerations
figure(); 
plot(x,aRP(x),x,aRP(x)-w^2*x);
grid on
title('accel by rad pressure (blue), total (red)')
%Clearly, the point of zero total accel is no longer at x=0, because of
%average thermal rad pressure! This could in principle be subtracted by
%shifting the spring equilibrium...

%% Simple ODE solver for a fixed initial condition at rest at x=0

%We start at x=0, p=0. If no rad pressure then no motion there. Any motion
%thus comes from rad pressure.
x0p0 = [0;0];
%However, part of it is because of shifted spring zero. So we compare to
%evolution under average rad pressure if both valves always open.

%In the simulation, we cannot just use instantaneous accel. as that would
%be like a conservative force. We must use explicit thermalization
%dynamics. Here we ignore the noise and just use a noiseless occupation
%that has reaction rate kappa
% dn/dt = -kappa* ( (fH+fC)*n + nH*fH + nC*fC )
% (Violates fluct-diss thm, of course, but here just qualitative behavior)
% (In practice, noise contribution can be huge, however! Need to check...)
kappa = 10; %should be > w for good operation

%y(1) is spring x, y(2) is dx/dt, y(3) is n(t) of work mode
fODE = @(t,y)( [y(2); r^2*y(3) - w^2*y(1) - gamma*y(2); ...
    kappa*( (nH-y(3)).*fH(y(1)) + (nC-y(3)).*fC(y(1)) )] );
%avg. rad press, no dynamical occupation
fODE0 = @(t,y)( [y(2); r^2*(nH+nC)/2 - w^2*y(1) - gamma*y(2)] );

tmax = 50*2*pi;

optt = odeset('RelTol',1e-8);

tic
[t,y] = ode113(fODE, [0 tmax], [x0p0;nC],optt);
[t0,y0] = ode113(fODE0, [0 tmax], x0p0,optt);
toc

figure(); 
plot(t/2/pi*w,y(:,1),t0/2/pi*w,y0(:,1));
grid on
title('Spring position (blue), only avg press. (red)')
xlabel('time in units 2\pi/\omega (spring period)')

%1st conclusion:
% gain is there and excites the spring oscillation. but nontrivial
% dependence on r and kappa, as expected.
%Also, there is a modulation of the eigenfrequency when compared to the
%unmodulated case. As in optomech. cavity cooling there should be sidebands
%in the Fourier spectrum that explain this...