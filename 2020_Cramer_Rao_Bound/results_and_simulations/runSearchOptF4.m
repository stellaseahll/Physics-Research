%clear;clc;
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
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        [x1, f1] = myMinSearch(@(theta) findFI5optBRUTEFORCE(theta,gammat(i),nbar(j)),-0.1+b1_F4theta(i,j),0.1+b1_F4theta(i,j),1e-8);
%         [x2, f2] = myMinSearch(@(theta) findFI5optBRUTEFORCE(theta,gammat(i),nbar(j)),pi/2,pi,1e-8);
        b1_F5theta(i,j) = x1;
        b1_F5(i,j) = -f1;
    end
    save('optAngleF5.mat');
end

load('b1FN.mat');
gammat = linspace(0,1,101);
nbar = linspace(0,10,101);
gammat(1) = 0.0001;
nbar(1) = 0.001;
for i= 1:length(gammat)
    for j = 1:length(nbar)
        fprintf('(%d,%d)\n',i,j);
        [x1, f1] = myMinSearch(@(theta) findFI6optBRUTEFORCE(theta,gammat(i),nbar(j)),-0.1+b1_F4theta(i,j),0.1+b1_F4theta(i,j),1e-6);
%         [x2, f2] = myMinSearch(@(theta) findFI5optBRUTEFORCE(theta,gammat(i),nbar(j)),pi/2,pi,1e-8);
        b1_F6theta(i,j) = x1;
        b1_F6(i,j) = -f1;
    end
    save('optAngleF6.mat');
end