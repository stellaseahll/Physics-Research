function y = findActualOptc(theta)
%finds optimum theta that maximises energy 
%assumes c=1 but need not be the case
incE = MERunClassical(0,theta,300,50);
y = (-incE(1,end)+incE(1,end-1))/theta;
% fprintf('(%.3f,%.3f)\n',q,y);
end

