function y = findActualOptcErg(theta)
%finds optimum theta that maximises ergotropy 
%assumes c=1 but need not be the case
incE = MERunClassical(0,theta,300,50);
y = (-incE(2,end)+incE(2,end-1))/theta;
% fprintf('(%.3f,%.3f)\n',q,y);
end

