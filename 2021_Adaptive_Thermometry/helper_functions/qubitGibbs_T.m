function T = qubitGibbs_T(p)
%function to compute T of given excitation prob p in 0..1 and degeneracy d
%number of excited levels

T = 1./log(d./p - d);

T(p==0) = 0;
T(p==0.5) = inf;
T(p==1) = 0; %careful, branch cut!