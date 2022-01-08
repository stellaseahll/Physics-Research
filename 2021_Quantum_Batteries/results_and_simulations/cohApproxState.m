function P = cohApproxState(n0,theta,q,c,alpha,n,k)

phi = linspace(-pi,0,10000);
omega = sqrt(q*(1-q))*sin(2*theta);
ptheta = sin(theta)^2;
nu = ptheta*(1-2*q);
gamma = nu-c*omega*sin(phi-alpha);
mu = n0 + k*gamma;
sigma = sqrt(k*ptheta-k*gamma.^2);
f = exp(-0.5*((n-mu)./sigma).^2)./sigma/(2*pi)^1.5;
P = trapz(phi,f);
