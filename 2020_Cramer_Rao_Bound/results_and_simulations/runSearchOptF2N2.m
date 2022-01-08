clf;clear;clc;
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
gammat = linspace(0,1,101);
nbar = logspace(-4,0,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;
r = 0.5;
phi = 0;
for i= 1:length(gammat)
     fprintf('(%d)\n',i);
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
%         F2gg(i,j) = -findFI2N2opt([0;0],r,phi,gammat(i),nbar(j));
%         F2gp(i,j) = -findFI2N2opt([0;pi/2],r,phi,gammat(i),nbar(j));
%         F2pg(i,j) = -findFI2N2opt([pi/2;0],r,phi,gammat(i),nbar(j));
        [x1, f1] = fminsearch(@(x) findFI2N2opt(x,r,phi,gammat(i),nbar(j)),[0.1; 1.6]);
%         if (f1<f2)
%             F2(i,j) = -f1;
%             theta1(i,j) = x1(1);
%             theta2(i,j) = x1(2);
%             r(i,j) = 1;
%         else
%             F2(i,j) = -f2;
%             theta1(i,j) = x2(1);
%             theta2(i,j) = x2(2);
%             r(i,j) = 0.5;
%         end
        theta1(i,j) = x1(1);
        theta2(i,j) = x1(2);
        F2(i,j) = -f1;
        
        
    end
	save('optAngleF2N2_logspace1.mat');
end
