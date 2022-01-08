q = linspace(0,0.5,51);
theta = linspace(0,pi/2,51);
theta(1) = [];
theta(end) = [];
N = 500;
n = (0:N).';
n0 = 0;
rho0 = zeros(N+1); rho0(n0+1,n0+1) = 1;
for i = 1:length(q)
    for j = [1 length(theta)]
        fprintf('(%d,%d)\n',i,j);
        omega = sqrt(q(i)*(1-q(i)))*sin(theta(j)*2);
        v = sin(theta(j))^2*(1-2*q(i));
        ndotanatmp(i,j) = omega + v;
        vtmp = [1 2];
        tmax = 200;
        nt = 0;
        L = MERunQuantumGetL(q(i),theta(j),N);
        rcoh = rho0(:);
        for k=1:10
            rcoh = rcoh + L*rcoh; 
            nt(k+1) = sum(n.*rcoh(1:N+2:end));
        end
        vtmp = diff(nt(end-10:end)); % 5 points
        
        while range(vtmp)/mean(vtmp) >10^-5
            rcoh = rcoh + L*rcoh; 
            nt(end+1) =  sum(n.*rcoh(1:N+2:end));
            vtmp = diff(nt(end-10:end));
        end
        ndottmp(i,j) = mean(vtmp);
    end
end
save('data_contourSpeed2.mat');
% %                   
% while 1
%         i = randi(51);
%         j = randi(49);
%         fprintf('q=%.3f,theta=%.3f\n',q(i),theta(j)/pi);
%         omega = sqrt(q(i)*(1-q(i)))*sin(theta(j)*2);
%         v = sin(theta(j))^2*(1-2*q(i));
%         ndotana(i,j) = omega + v;
%         vtmp = [1 2];
%         tmax = 200;
%         nt = 0;
%         L = MERunQuantumGetL(q(i),theta(j),N);
%         rcoh = rho0(:);
%         for k=1:10
%             rcoh = rcoh + L*rcoh; 
%             nt(k+1) = sum(n.*rcoh(1:N+2:end));
%         end
%         vtmp = diff(nt(end-10:end)); % 5 points
%         
%         while range(vtmp)/mean(vtmp) >10^-5
%             rcoh = rcoh + L*rcoh; 
%             nt(end+1) =  sum(n.*rcoh(1:N+2:end));
%             vtmp = diff(nt(end-10:end));
%         end
%         ndot(i,j) = mean(vtmp);
%         subplot(1,2,1);
%         plot(0:length(nt)-2,diff(nt),'k',[0 length(nt)-2],ndot(i,j)*[1 1],'bo-');
%         subplot(1,2,2);
%         plot(0:length(nt)-1,nt,'k',[0 length(nt)-1],ndot(i,j)*[0 length(nt)-1],'bo-');
%         pause();
% end
%         