function y = findOptq(q,n)
%given n number of probes over an interval pi/2, find the optimal q to maximise drift velocity
%assumes c=1 but need not be the case
theta = pi/2/n;
omega = sqrt(q*(1-q))*sin(2*theta);
v = sin(theta)^2*(1-2*q);
y = -n*(v+omega); %negative because maximising with fminsearch
end

