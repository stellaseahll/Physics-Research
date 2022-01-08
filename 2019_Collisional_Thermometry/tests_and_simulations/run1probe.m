nbar = 10*ones(1,N);
N = 10000;
gt = rand(1,N)*0.01*pi;
gammat = rand(1,N);
for i = 1:N
    fprintf('(%d,%d)\n',i);
    s = spSinglePassXY(gt(i),gammat(i),nbar(i),nbar(i)*0.0001,[1 0; 0 0],6,1);
    Fg(i,:) = s.getAllFish(); 
    s = spSinglePassXY(gt(i),gammat(i),nbar(i),nbar(i)*0.0001,[0 0; 0 1],6,1);
    Fe(i,:) = s.getAllFish(); 
end

save('randSmallG,nbar10.mat');

nbar = 1*ones(1,N);
N = 10000;
gt = rand(1,N)*0.01*pi;
gammat = rand(1,N);

for i = 1:N
    fprintf('(%d,%d)\n',i);
    s = spSinglePassXY(gt(i),gammat(i),nbar(i),nbar(i)*0.0001,[1 0; 0 0],6,1);
    Fg(i,:) = s.getAllFish(); 
    s = spSinglePassXY(gt(i),gammat(i),nbar(i),nbar(i)*0.0001,[0 0; 0 1],6,1);
    Fe(i,:) = s.getAllFish(); 
end

save('randSmallG,nbar1.mat');

nbar = 0.1*ones(1,N);
N = 10000;
gt = rand(1,N)*0.01*pi;
gammat = rand(1,N);

for i = 1:N
    fprintf('(%d,%d)\n',i);
    s = spSinglePassXY(gt(i),gammat(i),nbar(i),nbar(i)*0.0001,[1 0; 0 0],6,1);
    Fg(i,:) = s.getAllFish(); 
    s = spSinglePassXY(gt(i),gammat(i),nbar(i),nbar(i)*0.0001,[0 0; 0 1],6,1);
    Fe(i,:) = s.getAllFish(); 
end

save('randSmallG,nbar0.1.mat');
% nbar = 1.0;
% gt = 0.1:0.1:2;
% gt(gt==1) = [];
% gt(gt==2) = [];
% gt = gt*pi;
% gammat = logspace(-3,2,20);
% 
% for i = 1:length(gt)
%     for j = 1:length(gammat)
%         fprintf('(%d,%d)\n',i,j);
%         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 0; 0 0],6,1);
%         F = s.getAllFish();
%         Fg(i,j) = F(end);
%         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1+1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1-1/sqrt(2)]/2,2,1);
%         F = s.getAllFish(); 
%         F45(i,j) = F(end);
% %         s.alg = 3;
% %         F(j,:) = s.getAllFish(); hold on;
% %         Fg1(i,j) = k(1);
% %         Fg2(i,j) = k(end);
% %         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[0 0; 0 1],2,1);
% %         s.alg = 2;
% %         k = s.getAllFish();
% %         Fe1(i,j) = k(1);
% %         Fe2(i,j) = k(end);
% %         Fe(i,j) = k(end);
% %         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 1; 1 1]/2,2,1);
% %         s.alg = 3;
% %         k = s.getAllFish(); hold on;
% %         Fp1(i,j) = k(1);
% %         Fp2(i,j) = k(end);
% %         Fp(i,j) = k(end);
% %         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1+1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1-1/sqrt(2)]/2,2,1);
% %         s.alg = 3;
% %         k = s.getAllFish(); hold on;
% %         F45(i,j) = k(end);
%     end
% end
% 
% save('1probeNbar1.mat');
% 
% nbar = 10;
% gt = (0:0.01:1)*pi;
% gt([1 end]) = [];
% gammat = logspace(-3,2,100);
% 
% for i = 1:length(gt)
%     for j = 1:length(gammat)
%         fprintf('(%d,%d)\n',i,j);
%         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 0; 0 0],2,1);
%         s.alg = 3;
%         Fg(i,j) = s.getAllFish(); hold on;
%         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[0 0; 0 1],2,1);
%         s.alg = 3;
%         Fe(i,j) = s.getAllFish(); hold on;
%         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1 1; 1 1]/2,2,1);
%         s.alg = 3;
%         Fp(i,j) = s.getAllFish(); hold on;
%         s = spSinglePassXY(gt(i),gammat(j),nbar,nbar*0.0001,[1+1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1-1/sqrt(2)]/2,2,1);
%         s.alg = 3;
%         F45(i,j) = s.getAllFish(); hold on;
%     end
% end
% 
% save('1probeNbar10.mat');
% 
% nbar = 0.1;
% gt = (0:0.01:1)*pi;
% gt([1 end]) = [];
% gammat = logspace(-3,2,100);
% 
% 
