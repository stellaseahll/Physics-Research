%optimal classical parameters
[optCAngle optCdE] = fminsearch(@(theta) findActualOptcErg(theta),pi/2);
optCdE = -optCdE;

%actual simulation
N = 200; %battery level
numC = 50;
Theta = numC*optCAngle; %total angle Theta 
n = [2 4]; %number of probes used 
numQ = numC*n; %total number of quantum probes
theta = optCAngle./n; %angle per Q interaction
% q = 0:0.01:0.5;
% for j = 1:length(q)
%     for i = 1:length(n)
%         fprintf('(%d,%d)\n',i,j);
%         [cohE{i,j},incE{i,j},~,~] = MERun(q(j),theta(i),N,numQ(i));
%         lastcohE(i,j) = cohE{i,j}(end);
%         lastincE(i,j) = incE{i,j}(end);
%     end
% %     save('dataCompareClassical.mat');
% end

% [~, optqidx] = max(lastcohE');
qstart = 0;
qend = 0.5;
% q2(:,[1 end]) = [];
% clear cohE incE lastcohE lastincE
for i = 1:length(n)
    qopt(i) = myMinSearch(@(q) findActualOptqErg(q,theta(i),N,numQ(i)),qstart,qend,1e-6);
end
clear cohE incE lastcohE lastincE
for i = 1:length(n)
    fprintf('(%d,%d)\n',i,j);
    [cohE{i},incE{i},~,~] = MERun(qopt(i),theta(i),N,numQ(i));
    save('dataCompareClassicalErgotropy.mat');
end