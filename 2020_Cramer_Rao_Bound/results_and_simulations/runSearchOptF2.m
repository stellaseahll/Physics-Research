clear;clc;
load('b1FN.mat');
% gammat = linspace(0,1,501);
% gammat = gammat(2:end);
% nbar = linspace(0,10,501);
% nbar = nbar(2:end);
gammat = linspace(0,1,101);
nbar = linspace(0,10,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;

for i= 1:length(gammat)
     fprintf('(%d)\n',i);
    for j = 1:length(nbar)
%         fprintf('(%d,%d)\n',i,j);
        precision = 10^min([floor(log10(b1_F1theta(i,j)))-4,-8]);
        precision = max([1e-9 precision]);
            
        [x1, f1] = myMinSearch(@(theta) findFI2opt(theta,gammat(i),nbar(j)),0,pi/2,precision);
        f1 = -f1;
        if (f1>=2*b1_F1(i,j))
            b1_F2theta(i,j) = x1;
            b1_F2(i,j) = f1;
        else
            b1_F2theta(i,j) = b1_F1theta(i,j);
            b1_F2(i,j) = 2*b1_F1(i,j);
        end
    end
end
save('b1F2.mat');
load('b1FN.mat');
for i= 1:length(gammat)
%      fprintf('(%d)\n',i);
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        precision = 10^min([floor(log10(b1_F1(i,j)))-4,-9]);
        [x1, f1] = myMinSearch(@(theta) findFI4optBRUTEFORCE(theta,gammat(i),nbar(j)),b1_F4theta(i,j)-0.2,b1_F4theta(i,j)+0.2,precision);
        if (-f1>4*b1_F1(i,j))
            b1_F4theta(i,j) = x1;
            b1_F4(i,j) = -f1;
        else
            b1_F4theta(i,j) = b1_F1theta(i,j);
            b1_F4(i,j) = 4*b1_F1(i,j);
        end
    end
end
save('b1F4.mat');