function y = findActualOptq(q,theta,N,numQ)
%given n number of probes over an interval pi/2, find the optimal q to maximise drift velocity
%assumes c=1 but need not be the case

 [cohE,~,~,~] = MERun(q,theta,N,numQ);
 y = -cohE(1,end);
 fprintf('(%.3f,%.3f)\n',q,y)
end

